/* 
	Copyright 2012 Jeff Poole.  All rights reserved.
   	Portions of the code are close enough to the zfec code that the zfec license may apply (GNU GPL v2 and others) -- for more information, see:	https://tahoe-lafs.org/trac/zfec/browser/zfec/README.rst
   */

package jfec

import (
	"fmt"
	"os" // for Exit/Open
	"bytes" // for Compare
	"encoding/binary" // for binary.Write
	"runtime" // for parallel processing
	"runtime/debug" // for Stack / PrintStack
	"bufio" // for Reader
)

/*
	This code is a first stab at an implementation of the zfec library in Go.  
	It was inspired by the flŭd project, and was used as a way to write some "real" code in Go (even though much of it was straight translations of C code -- including comments, to help me keep my place), and as a way to remind myself how basic erasure codes work.
*/

////////// HELPER FUNCTIONS /////////////////////////////////

/* The smallest integer k such that k*d >= n. */
func div_ceil(n, d int) int {
	retval := (n/d)
	if n%d != 0 { retval = retval + 1 }
    return  retval
}

/*  The smallest integer k such that b^k >= n.

    log_ceil(n, 2) is the number of bits needed to store any of n values, e.g.
    the number of bits needed to store any of 128 possible values is 7. 
*/
func log_ceil(n, b uint) uint {
    p := uint(1)
    k := uint(0)
    for p < n {
        p *= b
        k += 1
    }
    return k
}


func assert(b bool) {
	if !b {
		fmt.Println("Assertion error!")
		debug.PrintStack()
		os.Exit(-1)
	}
}
///////////////////////////////////////////////////////////

type Fec struct {
	k, n uint8 /* in some places, n is called m... */
	enc_matrix []uint8 
}

func (f *Fec) String() string {
	return fmt.Sprintf("FEC: [k: %d n/m: %d]",f.k, f.n)
}


/*
 * Primitive polynomials - see Lin & Costello, Appendix A,
 * and  Lee & Messerschmitt, p. 453.
 */
const Pp = "101110001"
const STRIDE=8192
const CHUNKSIZE=4096


/*
 * To speed up computations, we have tables for logarithm, exponent and
 * inverse of a number.  We use a table for multiplication as well (it takes
 * 64K, no big deal even on a PDA, especially because it can be
 * pre-initialized an put into a ROM!), otherwhise we use a table of
 * logarithms. In any case the macro gf_mul(x,y) takes care of
 * multiplications.
 */

var gf_exp [510]uint8
var gf_log [256]int32
var inverse [256]uint8

/*
 * modnn(x) computes x % GF_SIZE, where GF_SIZE is 2**GF_BITS - 1,
 * without a slow divide.
 */
func modnn(x uint32) uint8 {
	for x >= 255 {
		x = x-255
		x = (x >> 8) + (x & 255)
	}
	return uint8(x)
}

/*
 * gf_mul(x,y) multiplies two numbers.  It is much faster to use a
 * multiplication table.
 *
 * USE_GF_MULC, GF_MULC0(c) and GF_ADDMULC(x) can be used when multiplying
 * many numbers by the same constant. In this case the first call sets the
 * constant, and others perform the multiplications.  A value related to the
 * multiplication is held in a local variable declared with USE_GF_MULC . See
 * usage in _addmul1().
 */
var gf_mul_table [256][256]uint8

func gf_mul(x,y uint8) uint8 {return gf_mul_table[x][y]}
var __gf_mulc_ [256]uint8 
func USE_GF_MULC() { return }
func GF_MULC0(c uint8) { __gf_mulc_ = gf_mul_table[c] }
func GF_ADDMULC(dst, x uint8)  { dst = dst ^ __gf_mulc_[x] }

/*
 * Generate GF(2**m) from the irreducible polynomial p(X) in p[0]..p[m]
 * Lookup tables:
 *     index->polynomial form		gf_exp[] contains j= \alpha^i;
 *     polynomial form -> index form	gf_log[ j = \alpha^i ] = i
 * \alpha=x is the primitive element of GF(2^m)
 *
 * For efficiency, gf_exp[] has size 2*GF_SIZE, so that a simple
 * multiplication of two numbers can be resolved without calling modnn
 */


func _init_mul_table() {
  	var i, j uint 
  	for i = 0; i < 256; i++ {
      	for j = 0; j < 256; j++ {
        	gf_mul_table[i][j] = gf_exp[modnn(uint32(gf_log[i] + gf_log[j]) )]
    	}
	}

  	for j = 0; j < 256; j++ {
      	gf_mul_table[0][j] = 0
      	gf_mul_table[j][0] = 0
  	}
}

func generate_gf() {
	var mask uint8 = 1
	gf_exp[8]=0 // I don't really understand this
    /*
     * first, generate the (polynomial representation of) powers of \alpha,
     * which are stored in gf_exp[i] = \alpha ** i .
     * At the same time build gf_log[gf_exp[i]] = i .
     * The first 8 powers are simply bits shifted to the left.
     */
	for i := 0; i < 8; i,mask = i+1,mask <<1 {
		gf_exp[i] = mask
		gf_log[gf_exp[i]] = int32(i)
        /*
         * If Pp[i] == 1 then \alpha ** i occurs in poly-repr
         * gf_exp[8] = \alpha ** 8
         */
		if Pp[i] == '1' {
			gf_exp[8] = gf_exp[8] ^ mask 
		}
	}
	/*
     * now gf_exp[8] = \alpha ** 8 is complete, so can also
     * compute its inverse.
     */
	gf_log[gf_exp[8]]=8
    /*
     * Poly-repr of \alpha ** (i+1) is given by poly-repr of
     * \alpha ** i shifted left one-bit and accounting for any
     * \alpha ** 8 term that may occur when poly-repr of
     * \alpha ** i is shifted.
     */
	mask = 1 << 7
	for i := 9; i < 255; i++ {
		if gf_exp[i-1] >= mask {
			gf_exp[i] = gf_exp[8] ^ ((gf_exp[i-1] ^ mask) << 1)
		} else {
			gf_exp[i] = gf_exp[i-1] << 1
		}
		gf_log[gf_exp[i]] = int32(i)
	}
	/*
     * log(0) is not defined, so use a special value
     */
	gf_log[0]=255 
	
	// set extended gf_exp values for fast multiply [gf_exp[n] == gf_exp[n%255]]
	for i :=0 ; i < 255 ; i++ {
		gf_exp[i+255] = gf_exp[i] 
	}
    /*
     * again special cases. 0 has no inverse. This used to
     * be initialized to 255, but it should make no difference
     * since noone is supposed to read from here.
     */
	inverse[0] = 0
	inverse[1] = 1
	for i := 2; i<= 255; i++ {
		inverse[i] = gf_exp[255-gf_log[i]]
	}
}
/*
 * Various linear algebra operations that i use often.
 */

 /*
 * addmul() computes dst[] = dst[] + c * src[]
 */
func addmul(dst, src []uint8, c uint8, sz uint) {
	if c == 0 {
		return // does nothing, so skip
	}
	var i uint 
	for i=0; i<sz; i = i + 1 {
		dst[i]=dst[i] ^ gf_mul_table[c][src[i]]
	}
}


/*
 * computes C = AB where A is n*k, B is k*m, C is n*m
 */

func _matmul(a, b, c []uint8, n,k,m uint8) {
	var row,col,i uint8
	var acc uint8
	for row = 0; row < n; row ++ {
		for col = 0; col < m; col++ {
			pa := row * k
			pb := col
			acc = 0
			for i=0; i<k; i,pa,pb = i+1,pa+1,pb+m {
				acc = acc ^ gf_mul(a[pa],b[pb])
			}
			c[row * m + col] = acc 
		}
	}
}

/*
 * _invert_mat() takes a matrix and produces its inverse
 * k is the size of the matrix.
 * (Gauss-Jordan, adapted from Numerical Recipes in C)
 * Return non-zero if singular.
 */
func _invert_mat(src []uint8, k uint) {
	var c uint8
	//var p *uint8
	var irow uint = 0
	var icol uint = 0
	var row,col,ix uint  // removed i because not used
	indxc :=  make([]uint,k)
	indxr :=  make([]uint,k)
	ipiv :=  make([]uint,k)
	id_row := make([]uint8,k) // "1*k" -- initialized to 0s
	/*
     * ipiv marks elements already used as pivots.
     */
    fmt.Printf("src matrix:\n% x\n",src)
    // C code zeros ipiv here, but that is automatic in Go
    for i:=0; i<int(k); i++ { ipiv[i]=0 }
    for col = 0; col < k; col++ {
    	var pivot_row []uint8 
        /*
         * Zeroing column 'col', look for a non-zero element.
         * First try on the diagonal, if it fails, look elsewhere.
         */
    	if ipiv[col] != 1 && src[col * k + col] != 0 {
    		irow = col
    		icol = col
    		fmt.Printf("Jumping from first check\n")
    		goto found_piv
    	}
    	for row = 0; row < k; row++ {
    		if ipiv[row] != 1 {
    			for ix = 0; ix < k; ix++ {
    				if ipiv[ix] == 0 {
    					if src[row*k+ix] != 0 {
    						irow = row
    						icol = ix
    						fmt.Printf("Jumping from second check\n")
    						goto found_piv
    					}
    				} else {
    					fmt.Printf("  ipiv[%d]=%d\n",ix,ipiv[ix])
    					assert(ipiv[ix] <= 1)
    				}
    			}
    		}
    	}
    	fmt.Printf("To found_piv by default?\n")
    found_piv:
    	fmt.Printf("  incrementing ipiv[%d] from %d\n",icol,ipiv[icol])
    	ipiv[icol] = ipiv[icol] + 1
    	/*
         * swap rows irow and icol, so afterwards the diagonal
         * element will be correct. Rarely done, not worth
         * optimizing.
         */
     	if irow != icol {
    	 	for ix=0; ix<k; ix++ {
         		src[irow*k + ix], src[icol*k + ix] = src[icol*k + ix], src[irow*k + ix] // SWAP(...)
         	}
        }
        indxr[col] = irow
        indxc[col] = icol
        pivot_row = src[icol*k:(icol+1)*k]
        c = pivot_row[icol]
        assert(c != 0)
        if c != 1 { // if c==1, this is a NOP
            /*
             * this is done often , but optimizing is not so
             * fruitful, at least in the obvious ways (unrolling)
             */
            c = inverse[c]
            pivot_row[icol] = 1
            for ix = 0; ix < k; ix++ {
            	pivot_row[ix] = gf_mul(c, pivot_row[ix])
            }
        }
        /*
        * from all rows, remove multiples of the selected row
        * to zero the relevant entry (in fact, the entry is not zero
        * because we know it must be zero).
        * (Here, if we know that the pivot_row is the identity,
        * we can optimize the addmul).
        */
		id_row[icol] = 1
		assert(len(pivot_row)==int(k) && len(id_row)==int(k))
		if bytes.Compare(pivot_row,id_row)!=0 {
			for ix=0; ix < k; ix++ {
				if ix != icol {
					c = src[(ix*k) + icol]
					src[(ix*k) + icol] = 0
					addmul(src[ix*k:(ix+1)*k], pivot_row, c, k)
				}
			}
		}
		id_row[icol] = 0
    }
    for col=k; col > 0; col-- {
    	if indxr[col-1] != indxc[col -1] {
    		for row=0; row < k; row++ {
    			src[row * k + indxr[col-1]], src[row * k + indxc[col-1]] = src[row * k + indxc[col-1]], src[row * k + indxr[col-1]]  // SWAP(...)
    		}
    	}
    }
}



/*
 * fast code for inverting a vandermonde matrix.
 *
 * NOTE: It assumes that the matrix is not singular and _IS_ a vandermonde
 * matrix. Only uses the second column of the matrix, containing the p_i's.
 *
 * Algorithm borrowed from "Numerical recipes in C" -- sec.2.8, but largely
 * revised for my purposes.
 * p = coefficients of the matrix (p_i)
 * q = values of the polynomial (known)
 */
func _invert_vdm(src []uint8, k uint8) {
	var i,j,row,col uint8

	//var b,c,p []uint8
	var t,xx uint8

	if k==1 {  /* degenerate case, matrix must be p^0 = 1 */
		return
	}
	/*
     * c holds the coefficient of P(x) = Prod (x - p_i), i=0..k-1
     * b holds the coefficient for the matrix inversion
     */
	c := make([]uint8,k)
	b := make([]uint8,k)
	p := make([]uint8,k)
	for j,i=1,0; i<k; i,j=i+1,j+k {
		c[i] = 0 // should be the default in Go
		p[i] = src[j]
	}
	/*
     * construct coeffs. recursively. We know c[k] = 1 (implicit)
     * and start P_0 = x - p_0, then at each stage multiply by
     * x - p_i generating P_i = x P_{i-1} - p_i P_{i-1}
     * After k steps we are done.
     */
	c[k-1] = p[0]
	for i = 1; i < k; i++ {
        p_i := p[i]            /* see above comment */
        for j = (k - 1 - (i - 1)); j < (k - 1); j++ {
            c[j] = c[j] ^ gf_mul(p_i, c[j + 1])
        }
        c[k - 1] = c[k - 1] ^ p_i
    }
    for row = 0; row < k; row++ {
        /*
         * synthetic division etc.
         */
        xx = p[row]
        t = 1
        b[k - 1] = 1 /* this is in fact c[k] */
        for i = (k - 1); i > 0; i-- {
            b[i-1] = c[i] ^ gf_mul(xx, b[i])
            t = gf_mul(xx, t) ^ b[i-1]
        }
        for col = 0; col < k; col++ {
            src[col * k + row] = gf_mul(inverse[t], b[col])
        }
    }
}

func init() {
	generate_gf()
	_init_mul_table()
}

/*
 * This section contains the proper FEC encoding/decoding routines.
 * The encoding matrix is computed starting with a Vandermonde matrix,
 * and then transforming it into a systematic matrix.
 */

const FEC_MAGIC = 0xFECC0DEC

func fec_free(f Fec) {  
	assert(false)  // no freeing in Go!
} 

// k is number of blocks required to reproduce data, m (n?) is the total number of blocks produced
func NewFec(k,n uint8) *Fec {
	var row,col uint8
	var i uint8
	assert(k<=n)
	// init {} always runs, so we don't have to check for init here
	retval := new(Fec)
	retval.k = k
	retval.n = n
	retval.enc_matrix = make([]uint8,n*k)
	tmp_m := make([]uint8,n*k)
	/*
     * fill the matrix with powers of field elements, starting from 0.
     * The first row is special, cannot be computed with exp. table.
     */
	tmp_m[0] = 1
	for col=1; col < k; col++ {
		tmp_m[col] = 0
	}
	for i,row=k,0; row < n-1; row,i = row+1,i+k {
		for col=0; col<k; col++ {
			tmp_m[i+col] = gf_exp[modnn(uint32(row)*uint32(col))]
		}
	}
	/*
     * quick code to build systematic matrix: invert the top
     * k*k vandermonde matrix, multiply right the bottom n-k rows
     * by the inverse, and construct the identity matrix at the top.
     */
    _invert_vdm (tmp_m, k)  /* much faster than _invert_mat */
    _matmul(tmp_m[k*k:len(tmp_m)], tmp_m[0:k*k], retval.enc_matrix[k*k:len(retval.enc_matrix)], n - k, k, k)
    /*
     * the upper matrix is I so do not bother with a slow multiply
     */
    for row = 0; row<k; row++ {
    	for col = 0; col < k; col++ {
    		if col==row {
    			retval.enc_matrix[row*k+col] = 1
    		} else {
    			retval.enc_matrix[row*k+col] = 0
    		}
    	}
    }
	return retval
}



/*
 * This generates a header that matches the zfec command-line tool.  The code is ugly, because it is a direct port.  Personally, I would spend the 2 extra bytes per file to make this code cleaner...
 */
func gen_header(code *Fec,padding, num uint) []byte {
	m := uint(code.n) // stupid n<->m
	k := uint(code.k)
	bitsused := uint(0)
	val := uint(0)
	val = val | (m-1)
	bitsused += 8 // first 8 bits aways encode m

	kbits := log_ceil(m,2) // number of bits needed to store all possible values of k
	val <<= kbits
	bitsused += kbits
	val = val | (k-1)

	padbits := log_ceil(k,2) // num bits needed to store all possible values of pad
	val <<= padbits
	bitsused += padbits
	val = val | padding

	shnumbits := log_ceil(m,2) //num bits needed to store all possible values of shnum
	val <<= shnumbits
	bitsused += shnumbits
	val = val | num 

	assert(bitsused >= 8)
	assert(bitsused <= 32)
	buf := new (bytes.Buffer)
	switch {
	case bitsused <= 16 :
		val <<= (16-bitsused)
		err := binary.Write(buf,binary.BigEndian, uint32(val))
		if err != nil {	fmt.Println("binary.Write failed:", err)  }
		return buf.Bytes()[2:4]
	case bitsused <= 24 :
		val <<= (24-bitsused)
		err := binary.Write(buf,binary.BigEndian, uint32(val))
		if err != nil {	fmt.Println("binary.Write failed:", err)  }
		return buf.Bytes()[1:4]
	}
	// default to 32-bit
	val <<= (32-bitsused)
	err := binary.Write(buf,binary.BigEndian, uint32(val))
	if err != nil {	fmt.Println("binary.Write failed:", err)  }
	return buf.Bytes()[0:4]
}



/* 
Helper function to take a filename and turn it into input and output bufio.Reader/Writer objects to pass to Encode_buffers() -- I wanted to make sure I could run encode on any arbitrary buffers, but the case of handling files makes it easier to compare to zfec. 
*/

func (code *Fec) Encode_files(baseFilename string, addZfecHeader bool) {
	infile,_ := os.Open(baseFilename)
	defer infile.Close()
	fi,_ := infile.Stat()
	insz := fi.Size()
	outfiles := make([]*bufio.Writer,code.n)
	padding := uint(uint(code.k)-uint(insz%int64(code.k)))
	for i:=0;i<int(code.n);i++ {
		fn := fmt.Sprintf("%s.%02d_%02d.jfec",baseFilename,i,code.n)
		fi,_ := os.Create(fn)
		defer fi.Close()
		outfiles[i] = bufio.NewWriter(fi)
		if addZfecHeader {
			outfiles[i].Write(gen_header(code,padding,uint(i)))
		}
	}
	code.Encode_buffers(bufio.NewReader(infile), outfiles, insz)
	// lets make sure everything is written to disk before we move on...
	for _,wr := range outfiles {
		wr.Flush()
	}
}

// Here I will try to use parallelism and sane IO (read some then write some) to make this work more smoothly with a smaller memory footprint
func (code *Fec) Encode_buffers(datain *bufio.Reader, dataout []*bufio.Writer, insz int64) {
	writerDone := make(chan bool,1)
	rateLimiter := make(chan int,runtime.NumCPU()+1) // no more than (number of CPUs+1) encoder processes running
	fmt.Printf("Setting GOMAXPROCS to NumCPU==%d\n",runtime.NumCPU())
	runtime.GOMAXPROCS(runtime.NumCPU()) // run on as many CPUs as I can

	kchunksz := int(CHUNKSIZE*int(code.k))
	numkchunks := int(insz / int64(kchunksz))
	if insz % int64(kchunksz) != 0 { numkchunks++ }

	completionChannels := make([]chan [][]uint8,numkchunks)
	for i,_ := range completionChannels {
		completionChannels[i] = make(chan [][]uint8,1)
	}
	curChan := 0

	go func() {
		for i,_ := range completionChannels {
			recdata := <-completionChannels[i]
			for j,wr := range dataout {
				wr.Write(recdata[j])
			}
		}

		writerDone<-true
	}()
	
	var blockLen,padding int 

	for  {
		// create buf every iteration so I can hand the blocks to goroutines and not worry about editing them after the fact
		buf := make([]byte,kchunksz) 
		numRead,_ := datain.Read(buf)
		if numRead == 0 {break} // jump out if done reading
		// create datablocks every iteration so I can hand the blocks to goroutines and not worry about editing them after the fact
		datablocks := make([][]uint8,code.k)
		if numRead == kchunksz {
			blockLen = CHUNKSIZE
			padding = 0
		} else {
			blockLen = div_ceil(numRead,int(code.k))
			padding = blockLen*int(code.k)-numRead
			// zero out padding...probably already zeroed, but the Go spec says it can use the extra space beyond the length it claims to have read
			for i:=numRead; i<len(buf) && i<numRead+padding; i++ { 	buf[i]=0 }
		}
		for i:=0; i<int(code.k); i++ {
			datablocks[i] = buf[i*blockLen:(i+1)*blockLen] //make([]uint8,blockLen)
		}
		rateLimiter <- 1
		go func(ch chan [][]uint8,chnum int) {
			fecs,dats := code.Fec_encode(datablocks,chnum)
			dats = append(dats,fecs...)
			ch<-dats
			<-rateLimiter
		}(completionChannels[curChan],curChan)
		curChan++
		
	}
	fmt.Println("Waiting for writer to get all results...")
	<-writerDone
	fmt.Printf("Done waiting! len(rateLimiter)=%d\n",len(rateLimiter))
}

func (code *Fec) Fec_encode(datablocks [][]uint8, chnum int) ([][]uint8,[][]uint8){
	var i,j uint
	var k uint 
	var stride uint 
	var blocksize uint

	blocksize = uint(len(datablocks[0]))

	block_nums := make([]uint8,code.n-code.k)
	
	fecs := make([][]uint8,code.n-code.k)
	for i=0;i<uint(code.n-code.k);i++ {
		block_nums[i]=code.k+uint8(i)
		fecs[i] = make([]uint8,blocksize)
	}
	//fmt.Printf("blocksize=%d, code.k=%d, code.n=%d\n\n",blocksize, code.k, code.n)
	for k=0; k<blocksize; k = k+STRIDE {
		if (blocksize-k) < STRIDE {
			stride = blocksize-k
		} else {
			stride = STRIDE 
		}
		//fmt.Printf("K=%d, stride=%d\n",k,stride)
		for i=0; int(i)<len(block_nums); i++ {
			fecnum := block_nums[i]
			assert(fecnum >= code.k)
			// fecs[i][k:k+stride] already zereoed by default
			for j=0; j<uint(code.k); j++ {
				addmul(fecs[i][k:k+stride], datablocks[j][k:k+stride], code.enc_matrix[uint(fecnum) * uint(code.k) + j], stride)
			}
		}
	}
	return fecs,datablocks
}


/**
 * Build decode matrix into some memory space.
 *
 * @param matrix a space allocated for a k by k matrix
 */
func build_decode_matrix_into_space(code *Fec, index []uint8,  k uint, matrix []uint8) {
    var i uint
    for i=0; i<k; i++ {
    	if uint(index[i]) < k {
    		// zero matrix[i*k:i*k+k], but that should be default...
    		matrix[i*k+i] = 1
    	} else {
    		copy(matrix[i*k:(i+1)*k], code.enc_matrix[uint(index[i]) * uint(code.k): uint(index[i]) * uint(code.k)+k])
    	}
    }
    _invert_mat(matrix, k)
}

/*
	I verified that the algorithm works with some basic tests, then left this as is.  It doesn't support zfec headers, doesn't run in parallel, etc. 
	Since the zfec source isn't very clear, inpkts is a slice of k uint8 arrays.  All the original data (block number n, where 0 <= n < k) need to be in their correct slot, if available.  Wherever original data is missing, you can subsitute one of the fec blocks (any block n, where k <= n < m).  The index slice holds the block numbers, so if k=3 and you have blocks 0 and 2, but are missing 1 which you replace with block 5, you might have index={0,5,2} with those block numbers matching up with the slices in inpkts.  The outpkts slice-of-slices only contains any rebuilt slices (in the case above, it will have one slice for n=1).  The sz parameter is the length of the slices you are passing in (and getting out).
*/
func (code *Fec) Fec_decode(inpkts, outpkts [][]uint8, index []uint8, sz uint) {
	m_dec := make([]uint8, code.k*code.k)
	var outix,row,col uint8
	outix = 0
	build_decode_matrix_into_space(code, index, uint(code.k), m_dec)
    for row=0; row<code.k; row++ {
        assert ((index[row] >= code.k) || (index[row] == row)) /* If the block whose number is i is present, then it is required to be in the i'th element. */
        if (index[row] >= code.k) {
            outpkts[outix] = make([]uint8, sz) // I think this will effectively zero it...if we need to
            for col=0; col < code.k; col++ {
                addmul(outpkts[outix], inpkts[col], m_dec[row * code.k + col], sz);
            }
            outix++
        }
    }	
}
package main

import (
	"jfec"
	"fmt"
	"os"
	"flag"
	)




func main() {
	var m,k int
	var addZfecHeader bool
	flag.IntVar(&m,"m",10,"total number of files to create")
	flag.IntVar(&k,"k",5,"minimum number of files needed to re-create original file")
	flag.BoolVar(&addZfecHeader,"zfecheader",true,"true if you want to add a header compatible with the zfec command line tool")
	flag.Parse()
	filenames := flag.Args()
	if len(filenames) == 0 {
		fmt.Println("At least one filename required on the command line to encode.")
		os.Exit(-1)
	}
	jf := jfec.NewFec(uint8(k),uint8(m))
	for _,fn := range filenames {
		fmt.Printf("Encoding file %s...\n",fn)
		jf.Encode_files(fn,addZfecHeader)
	}
}
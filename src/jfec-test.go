package main

import (
	"jfec"
	"fmt"
	"os"
	"flag"
	"runtime/debug"
	)


func assert(b bool) {
	if !b {
		fmt.Println("Assertion error!")
		debug.PrintStack()
		os.Exit(-1)
	}
}


func main() {
	var m,k int
	var addZfecHeader bool
	flag.IntVar(&m,"m",10,"total number of files to create")
	flag.IntVar(&k,"k",5,"minimum number of files needed to re-create original file")
	flag.BoolVar(&addZfecHeader,"zfecheader",true,"true if you want to add a header compatible with the zfec command line tool")
	flag.Parse()
	filenames := flag.Args()
	if len(filenames) 
	jf := jfec.NewFec(uint8(k),uint8(m))
	for fn := range filenames {
		fmt.Printf("Encoding file %s...\n",fn)
		jf.Encode_files(fn,addZfecHeader)
	}
}
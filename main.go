package main

import (
	"fmt"
	"strconv"
	"time"
)

// Program begins here
func main() {
	start := time.Now()
	// load map relating charmm to amoeba types for atoms in our bilayer
	fmt.Println("Reading CHARMM to AMOEBA atom type map...")
	convPath := "C:\\Users\\jtgou\\OneDrive\\Documents\\UT_Austin\\ren_lab\\lipids\\pren\\converter\\dmpg_type_conversion.txt"
	charmm2amoeba := retrieveCHARMM2AMOEBAmap(convPath)

	// batch convert all CHARMM bilayers in the in directory...
	inDir := "C:\\Users\\jtgou\\OneDrive\\Documents\\UT_Austin\\ren_lab\\lipids\\pren\\bilayersIn"
	// to AMOEBA bilayers in the out directory
	outDir := "C:\\Users\\jtgou\\OneDrive\\Documents\\UT_Austin\\ren_lab\\lipids\\pren\\bilayersOut"
	fmt.Println("Batch converting bilayers...")

	charmm2amoebaBatch(inDir, outDir, charmm2amoeba)
	duration := time.Since(start)
	fmt.Println("All bilayers converted.")
	fmt.Println("Total Runtime = " + strconv.Itoa(int(duration.Milliseconds())) + " ms.")

}






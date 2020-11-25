package main

import (
	"fmt"
	"log"
	"os"
	filepath2 "path/filepath"
	"strconv"
)

func writeFragment(atoms map[int]*atom, fragPath string) {
	err := os.MkdirAll(filepath2.Dir(fragPath), 0755)
	if err != nil {
		fmt.Println("Failed to create new directory: " + filepath2.Dir(fragPath))
		log.Fatal(err)
	}

	thisFile, err := os.Create(fragPath)
	if err != nil {
		fmt.Println("Failed to create new fragment file: " + fragPath)
		log.Fatal(err)
	}

	// write header
	_, err = thisFile.WriteString(strconv.Itoa(len(atoms)) + "\t" + filepath2.Base(fragPath) + "\n")
	if err != nil {
		fmt.Println("Failed to write header line to key: " + fragPath)
		log.Fatal(err)
	}

	// write body
	for i := 1; i <= len(atoms); i++ {

		line := strconv.Itoa(i) + "\t" + atoms[i].element + "\t" + fmt.Sprintf("%.6f", atoms[i].pos[0]) + "\t" +
				fmt.Sprintf("%.6f", atoms[i].pos[1]) + "\t" + fmt.Sprintf("%.6f", atoms[i].pos[2]) + "\t" +
				strconv.Itoa(atoms[i].atomType)

		for _, bondedAtom := range atoms[i].bondedAtoms {
			line += "\t" + strconv.Itoa(bondedAtom)
		}

		_, err = thisFile.WriteString(line + "\n")
		if err != nil {
			fmt.Println("Failed to write header line to key: " + fragPath)
			log.Fatal(err)
		}
	}
}
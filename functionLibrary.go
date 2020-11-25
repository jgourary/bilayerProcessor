package main

import (
	"bufio"
	"errors"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	filepath2 "path/filepath"
	"strconv"
	"strings"
	"sync"
)

// load map relating charmm to amoeba types for atoms in our bilayer
func retrieveCHARMM2AMOEBAmap(path string) map[string]string {
	// open file
	thisFile, err := os.Open(path)
	if err != nil {
		fmt.Println("Failed to open CHARMM to AMOEBA type conversion file: " + path)
		log.Fatal(err)
	}

	// Initialize scanner
	scanner := bufio.NewScanner(thisFile)

	charmm2amoeba := make(map[string]string)

	// iterate over all other lines
	for scanner.Scan() {
		// get next line
		line := scanner.Text()
		// split by whitespace
		tokens := strings.Fields(line)
		if len(tokens) > 1 {
			charmm2amoeba[tokens[0]] = tokens[1]
		}
	}
	return charmm2amoeba
}

// batch convert all CHARMM bilayers to AMOEBA bilayers in parallel
func charmm2amoebaBatch(inDir string, outDir string, charmm2amoeba map[string]string) {
	fileInfo, err := ioutil.ReadDir(inDir)
	if err != nil {
		fmt.Println("failed to read directory: " + inDir)
		log.Fatal(err)
	}
	err = os.MkdirAll(outDir, 0755)
	if err != nil {
		fmt.Println("failed to create dir: " + outDir)
		log.Fatal(err)
	}
	wg := sync.WaitGroup{}
	for i := 0; i < len(fileInfo); i++ {
		if filepath2.Ext(fileInfo[i].Name()) == ".pdb" {
			baseName := strings.Split(fileInfo[i].Name(),".")[0]
			pdbPath := filepath2.Join(inDir, fileInfo[i].Name())
			psfPath := filepath2.Join(inDir, baseName + ".psf")
			outPath := filepath2.Join(outDir, baseName + ".txyz")
			wg.Add(1)
			go charmm2intermediateSingle(pdbPath, psfPath, outPath, charmm2amoeba, &wg)
		}
	}
	wg.Wait()
}

// convert a single CHARMM bilayer to an AMOEBA bilayer
func charmm2intermediateSingle(pdb string, psf string, outPath string, charmm2amoeba map[string]string, wg *sync.WaitGroup) {
	// get atom information from pdb
	atoms := processPDB(pdb)
	// get bond information from psf
	bonds := processPSF(psf)

	// applies bond information from PSF to atom information from PDB
	for k, v := range bonds {
		atoms[k].bondedAtoms = v
	}

	postProcessBonds(atoms)
	postProcessAtomTypes(atoms)

	assignAMOEBAtypes(atoms, charmm2amoeba)

	writeFragment(atoms, outPath)
	wg.Done()
}

// assign correct AMOEBA types to all atoms from their charmm types
func assignAMOEBAtypes(atoms map[int]*atom, charmm2amoeba map[string]string) {
	for _, atom1 := range atoms {
		atomType, err := strconv.Atoi(charmm2amoeba[atom1.charmmID])
		if err != nil {
			log.Fatal("failed to convert atom type " + charmm2amoeba[atom1.charmmID] + " for charmm atom ID" + atom1.charmmID + " to integer")
		}
		atom1.atomType = atomType
	}
}

// corrects hydrogens bound to each other in psf file for water molecules in CHARMM bilayers
func postProcessBonds(atoms map[int]*atom) {
	for atomNum, atom1 := range atoms {
		for _, bondedAtomNum := range atom1.bondedAtoms {
			if atomNum < bondedAtomNum && atom1.element == "H" && atoms[bondedAtomNum].element == "H" {
				disconnect(atoms, atomNum, bondedAtomNum)
			}
		}
	}
}

// Corrects ion atom types with relevant +/-
func postProcessAtomTypes(atoms map[int]*atom) {
	for _, atom1 := range atoms {
		if atom1.element == "POT" { atom1.element = "K" }
	}
	for _, atom1 := range atoms {
		if len(atom1.bondedAtoms) == 0 {
			if atom1.element == "K" || atom1.element == "Na" || atom1.element == "Li" || atom1.element == "H" {
				atom1.element = atom1.element + "+"
			} else if atom1.element == "Cl" || atom1.element == "F" || atom1.element == "Br" || atom1.element == "I" {
				atom1.element = atom1.element + "-"
			}
		}
	}


}

// removes each atom from the other's list of bonded atoms
func disconnect(atoms map[int]*atom, atom1 int, atom2 int) {
	for i, bondedAtom := range atoms[atom1].bondedAtoms {
		if bondedAtom == atom2 {
			atoms[atom1].bondedAtoms = remove(atoms[atom1].bondedAtoms,i)
		}
	}
	for i, bondedAtom := range atoms[atom2].bondedAtoms {
		if bondedAtom == atom1 {
			atoms[atom2].bondedAtoms = remove(atoms[atom2].bondedAtoms,i)
		}
	}
}

// reads in CHARMM PDB and gets atom information
func processPDB(file string) map[int]*atom {

	// Create structure to store atoms
	atoms := make(map[int]*atom)

	// open file
	thisFile, err := os.Open(file)
	if err != nil {
		fmt.Println("Failed to open pdb file: " + file)
		log.Fatal(err)
	}

	// Initialize scanner
	scanner := bufio.NewScanner(thisFile)
	// iterate over all other lines
	for scanner.Scan() {
		// get next line
		line := scanner.Text()
		// split by whitespace
		tokens := strings.Fields(line)

		if len(tokens) > 10 {
			if tokens[0] == "ATOM" {



				// create new atom
				var newAtom atom

				// get number of atom from file
				atomNum, err := strconv.Atoi(tokens[1])
				if err != nil {
					newErr := errors.New("failed to convert token to an integer")
					log.Fatal(newErr)
				}

				// assign charmmID & element
				newAtom.charmmID = tokens[2]
				if tokens[2] == "POT" {
					newAtom.element = "K"
				} else {
					newAtom.element = string(tokens[2][0])
				}

				newAtom.molType = tokens[3]
				molNum, err := strconv.Atoi(tokens[4])
				newAtom.molNum = molNum

				// assign positions
				pos := make([]float64,3)
				for j := 5; j < 8; j++ {
					pos[j-5], err = strconv.ParseFloat(tokens[j],64)
					if err != nil {
						newErr := errors.New("Failed to convert token in position 0 on line " + strconv.Itoa(j) + " to a float64")
						log.Fatal(newErr)
					}
				}
				newAtom.pos = pos

				atoms[atomNum] = &newAtom

			}
		}
	}

	return atoms
}

// reads in CHARMM PSF and gets bond information
func processPSF(file string)  map[int][]int {
	// open file
	thisFile, err := os.Open(file)
	if err != nil {
		fmt.Println("Failed to open psf file: " + file)
		log.Fatal(err)
	}

	// Initialize scanner
	scanner := bufio.NewScanner(thisFile)

	isInBondsSection := false
	bonds := make(map[int][]int)

	// iterate over all other lines
	for scanner.Scan() {
		// get next line
		line := scanner.Text()
		// split by whitespace
		tokens := strings.Fields(line)

		if isInBondsSection {
			for i := 0; i < len(tokens)-1; i++ {
				if i % 2 == 0 {
					atom1, err := strconv.Atoi(tokens[i])
					atom2, err := strconv.Atoi(tokens[i+1])
					if err == nil {
						if !isInSlice(bonds[atom1], atom2) {
							bonds[atom1] = append(bonds[atom1], atom2)
						}
						if !isInSlice(bonds[atom2], atom1) {
							bonds[atom2] = append(bonds[atom2], atom1)
						}
					}

				}
			}
		}

		if len(tokens) == 3 {
			if tokens[1] == "!NBOND:" {
				isInBondsSection = true
			} else {
				isInBondsSection = false
			}
		}


	}
	return bonds
}

// removes element from slice
func remove(s []int, i int) []int {
	s[i] = s[len(s)-1]
	return s[:len(s)-1]
}

// checks if element is in slice
func isInSlice(slice []int, val int) bool {
	isPres := false
	for _, integer := range slice {
		if integer == val {
			isPres = true
		}
	}
	return isPres
}





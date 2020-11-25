package main

// Atom data type
type atom struct {
	element string
	atomType int
	bondedAtoms []int

	pos []float64

	charmmID string
	molNum int
	molType string
}

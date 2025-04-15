unit STDMol;

{$mode objfpc}{$H+}

interface

const
  natoms = 4096;

type
	TAtomRecord = record
	id: integer;
	atom_name: string [3];
	atom_type: string;
	charge: real;
	x: real;
	y: real;
	z: real;
end;

type
	TSTDMol = array [1..natoms] of TAtomRecord;

implementation

End.


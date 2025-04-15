unit Mol2;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, STDMol;

type
  TMol2AtomRecord = record
  id: integer;
	atom_id: ansistring;
	atom_name: string[3];
	x: real;
	y: real;
	z: real;
	atom_type: string;
	subst_id: string;
	subst_name: string;
	charge: real;
	stats_bit: string[2];
end;

type
  TTripos = record
  atoms: array[0..natoms] of TMol2AtomRecord;
  nat: integer;
end;

function mol2Process (lines: array of string): TTripos;

implementation

var
  mol: TTripos;

procedure process_atom (line: string; id: integer);
var
  i: integer = 0;
  col: integer = 0;
  chrs: ansistring;
  charge: string = '';
  x: string = '';
  y: string = '';
  z: string = '';
	setcol: boolean = false;

begin
  mol.atoms[id].id:= id;
  mol.atoms[id].atom_id:= '';
  mol.atoms[id].atom_name:= '';
  mol.atoms[id].atom_type:= '';
  mol.atoms[id].charge:= 0.0;
  mol.atoms[id].stats_bit:= '';
  mol.atoms[id].subst_id:= '';
  mol.atoms[id].x:= 0;
  mol.atoms[id].y:= 0;
  mol.atoms[id].z:= 0;
  mol.atoms[id].subst_name:= '';
  for i:= 0 to length(line) do begin
      chrs:= copy(line, i, 1);
      if not (chrs = chr(13)) then begin
          if ((chrs = ' ') or (chrs = chr(9))) then begin
              setcol:= false;
          end;
          if (not (chrs = ' ') and (setcol = false)) then begin
              setcol:= true;
              col:= col + 1;
          end;
          if ((setcol = true) and (col = 1)) then begin
          	mol.atoms[id].atom_id:= mol.atoms[id].atom_id + chrs;
          end;
          if ((setcol = true) and (col = 2)) then begin
          	mol.atoms[id].atom_name:= mol.atoms[id].atom_name + chrs;
          end;
          if ((setcol = true) and (col = 3)) then begin
            x:= x + chrs;
          end;
          if ((setcol = true) and (col = 4)) then begin
            y:= y + chrs;
          end;
          if ((setcol = true) and (col = 5)) then begin
            z:= z + chrs;
          end;
          if ((setcol = true) and (col = 6)) then begin
          	mol.atoms[id].atom_type:= mol.atoms[id].atom_type + chrs;
          end;
          if ((setcol = true) and (col = 7)) then begin
          	mol.atoms[id].subst_id:= mol.atoms[id].subst_id + chrs;
          end;
          if ((setcol = true) and (col = 8)) then begin
          	mol.atoms[id].subst_name:= mol.atoms[id].subst_name + chrs;
          end;
          if ((setcol = true) and (col = 9)) then begin
            charge:= charge + chrs;
          end;
          if ((setcol = true) and (col = 10)) then begin
          	mol.atoms[id].stats_bit:= mol.atoms[id].stats_bit + chrs;
          end;
      end;
  end;
  mol.atoms[id].x:= strtofloat(x);
  mol.atoms[id].y:= strtofloat(y);
  mol.atoms[id].z:= strtofloat(z);
  mol.atoms[id].charge:= strtofloat(charge);
end;

function mol2Process (lines: array of string): TTripos;
var
	line: string = '';
  rec: string = '';
  setrec: boolean;
  recn: integer = 0;
  id: integer = 0;
begin
  for line in lines do begin
    setrec:= false;
    if (copy(line, 0, 9) = '@<TRIPOS>') then begin
    	setrec:= true;
    end;
    if (setrec) then begin
      rec:= copy (line, 10, 12);
    end;
    if (setrec and (length(rec) > 0)) then begin
        recn:= recn + 1;
        continue;
    end;
    case rec of
    	'ATOM': begin process_atom (line, id); id:= id + 1; end;
    end;
  end;
  if (id = 0) then begin raise exception.create('ERROR in mol2 file'); end;
  mol.nat:= id;
  mol2Process:= mol;
end;



end.


unit Coord;

{$mode objfpc}{$H+}

interface

uses
Classes, SysUtils, TypeDef;

type
  TCoordinates = record
    atom: array of TVector;
    charge: array of real;
    weight: TVector;
    wname: string;
    atype: array of string;
    coeff: real;
  end;

implementation

end.


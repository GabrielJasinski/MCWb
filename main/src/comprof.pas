unit ComProF;

{$mode objfpc}{$H+}

interface

uses
Classes, SysUtils, Math, Typ, TypeDef;

function roundd (number: real; decimal: integer): extended;
function split (str: string; mark: char): TSVector;

implementation

function roundd (number: real; decimal: integer): extended;
var
	places: double;
begin
  places:= power(10, decimal);
  roundd:= round(number * places) / places;
end;

function split (str: string; mark: char): TSVector;
var
  i: integer = 0;
  j: integer = 0;
  chr: string = '';
  field: string = '';
  splits: TSVector;
begin
  for i:= 1 to length(str) do begin
  	chr:= copy(str, i, 1);
    if (chr = mark) then begin
    	setlength (splits, length(splits) + 1);
      splits[j]:= field;
      field:= '';
      j:= j + 1;
    end else begin
      field:= field + chr;
    end;
  end;
  split:= splits;
end;

end.


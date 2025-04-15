unit Weights;

{$mode objfpc}{$H+}

interface

uses Classes, SysUtils, StrUtils, TypeDef, OFileTextArray;

type
  TWeightsxl = array of TVector;
  TWeightList = array of string;
  TWcompiled = array of TVector;

  TWeights = record
    nwe: integer;
    nat: integer;
    weight: array of string;
    atype: array of string;
    molpar: array of boolean;
    value: array of array of real;
  end;

procedure loadWeights (f: string);
function getWeights(): TWeightList;
function getWeightvalue(at: string; w: integer): real;

implementation
var
  aweights: TWeights;
  fname: string = '';
  wlines: array of string;
  atypes: array of string;
  weightname: array of string;
  weighttype: array of string;
  weightvalue: array of array of real;
procedure loadWeights (f: string);
var
  i, x: integer;
  valstr: string;
begin
  fname:= f;
	wlines:= readFileText (fname);
 	aweights.nwe:= wordCount(wlines[1], [',']) - 1;
  aweights.nat:= length(wlines) - 2;
 	setLength (aweights.atype, aweights.nat);
  setLength (aweights.molpar, aweights.nat);
  setLength (aweights.weight, aweights.nwe);
  setLength (aweights.value, aweights.nat, aweights.nwe);
  for x:= 1 to wordCount(wlines[1], [',']) - 1 do begin
    aweights.weight[x-1]:= trim(extractDelimited(x + 1, wlines[1], [',']));
  end;
 	for i:= 0 to (length(wlines) - 1) do
 	begin
    for x:= 1 to (wordCount(wlines[i], [','])) do begin
  		if ((i > 1) and (x = 1)) then	begin
    		aweights.atype[i-2]:= trim(extractDelimited(x, wlines[i], [',']));
    	end;
    	if ((i > 1) and (x > 1)) then	begin
        valstr:= trim(extractDelimited(x, wlines[i], [',']));
        if ((valstr = 'coor') or (valstr = 'qnor') or (valstr = 'qstd')) then begin
          aweights.value[i-2][x-2]:= 0;
          aweights.molpar[i-2]:= true;
        end
        else
        begin
          aweights.value[i-2][x-2]:= strtofloat(valstr);
          aweights.molpar[i-2]:= false;
        end;
    	end;
    end;
  end;
end;
function getWeights(): TWeightList;
begin
  getWeights:= aweights.weight;
end;
function getWeightvalue(at: string; w: integer): real;
var
  i: integer;
  atfound: boolean = false;
begin
  for i:= 0 to length (aweights.atype) - 1 do begin
    if (aweights.atype[i] = at) then begin
      atfound:= true;
    	getWeightvalue:= aweights.value[i][w];
      exit;
    end;
  end;
  if atfound = false then begin
    raise exception.create(('** WARNING ** Atom Type [' + at + '] do not have parameters. This molecule is skipped. '));
  end;
  getWeightvalue:= 0;
end;
end.


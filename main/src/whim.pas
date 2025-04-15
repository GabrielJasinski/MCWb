unit WHIM;

{$mode objfpc}{$H+}

interface

uses
Classes, SysUtils, Math, Typ;

type
  TWHIMresult = record
  	descriptor: string;
    value: array of real;
  end;
  TEigenvalues = array[1..3] of complex;
  TEigenvectors =  array [1..3, 1..3] of complex;
  TWHIMEigenvalues = array[1..3] of real;
  TWHIMEigenvectors = array [1..3, 1..3] of ArbFloat;
  TPCAScores = array of array of real;


function whimLkw (eigval: TWHIMEigenvalues): TWHIMresult;
function whimPkw (eigval: TWHIMEigenvalues): TWHIMresult;
function whimTw (eigval: TWHIMEigenvalues): TWHIMresult;
function whimAw (eigval: TWHIMEigenvalues): TWHIMresult;
function whimVw (eigval: TWHIMEigenvalues): TWHIMresult;
function whimKw (eigval: TWHIMEigenvalues): TWHIMresult;
function whimGkw (pcascores: TPCAScores; nat: integer; symtol: real): TWHIMresult;
function whimGw (pcascores: TPCAScores; nat: integer; symtol: real): TWHIMresult;
function whimEkw (eigval: TWHIMEigenvalues; sqs: array of real; nat: integer): TWHIMresult;
function whimDw (eigval: TWHIMEigenvalues; sqs: array of real; nat: integer): TWHIMresult;

implementation

function whimLkw (eigval: TWHIMEigenvalues): TWHIMresult;
var
  whimres: TWHIMresult;
begin
  whimres.descriptor:= 'L';
  setlength (whimres.value, 3);
  whimres.value[0]:= eigval[1];
  whimres.value[1]:= eigval[2];
  whimres.value[2]:= eigval[3];
  whimLkw:= whimres;
end;

function whimPkw (eigval: TWHIMEigenvalues): TWHIMresult;
var
  whimres: TWHIMresult;
  eigvals: real;
begin
  whimres.descriptor:= 'P';
  setlength (whimres.value, 2);
  eigvals:= eigval[1] + eigval[2] + eigval[3];
  whimres.value[0]:= eigval[1] / eigvals;
  whimres.value[1]:= eigval[2] / eigvals;
  //whimres.value[2]:= eigval[3] / eigvals;
  whimPkw:= whimres;
end;

function whimEkw (eigval: TWHIMEigenvalues; sqs: array of real; nat: integer): TWHIMresult;
var
  whimres: TWHIMresult;
begin
  whimres.descriptor:= 'E';
  setlength (whimres.value, 3);
  whimres.value[0]:= 0;
  whimres.value[1]:= 0;
  whimres.value[2]:= 0;
  if (sqs[0] > 0) then begin whimres.value[0]:= (sqr(eigval[1]) * nat) / sqs[0]; end;
  if (sqs[1] > 0) then begin whimres.value[1]:= (sqr(eigval[2]) * nat) / sqs[1]; end;
  if (sqs[2] > 0) then begin whimres.value[2]:= (sqr(eigval[3]) * nat) / sqs[2]; end;
  whimEkw:= whimres;
end;

function whimTw (eigval: TWHIMEigenvalues): TWHIMresult;
var
  whimres: TWHIMresult;
begin
  whimres.descriptor:= 'T';
  setlength (whimres.value, 1);
  whimres.value[0]:= eigval[1] + eigval[2] + eigval[3];
  whimTw:= whimres;
end;

function whimAw (eigval: TWHIMEigenvalues): TWHIMresult;
var
  whimres: TWHIMresult;
begin
  whimres.descriptor:= 'A';
  setlength (whimres.value, 1);
  whimres.value[0]:= eigval[1]*eigval[2] + eigval[1]*eigval[3] + eigval[2]*eigval[3];
  whimAw:= whimres;
end;

function whimVw (eigval: TWHIMEigenvalues): TWHIMresult;
var
  whimres: TWHIMresult;
  eigvals: real;
begin
  whimres.descriptor:= 'V';
  setlength (whimres.value, 1);
  whimres.value[0]:= eigval[1] + eigval[2] + eigval[3] + (eigval[1]*eigval[2] + eigval[1]*eigval[3] + eigval[2] * eigval[3]) + eigval[1]*eigval[2]*eigval[3];
  whimVw:= whimres;
end;

function whimKw (eigval: TWHIMEigenvalues): TWHIMresult;
var
  whimres: TWHIMresult;
  eigvals: real;
  T: real;
begin
  whimres.descriptor:= 'K';
  setlength (whimres.value, 1);
  eigvals:= eigval[1] + eigval[2] + eigval[3];
  T:= (3/4) * (sqrt(sqr(((eigval[1] / eigvals) - (1/3)))) +  sqrt(sqr(((eigval[2] / eigvals) - (1/3)))) +sqrt(sqr(((eigval[3] / eigvals) - (1/3)))));
  whimres.value[0]:= T;
  whimKw:= whimres;
end;

function whimDw (eigval: TWHIMEigenvalues; sqs: array of real; nat: integer): TWHIMresult;
var
  whimres: TWHIMresult;
  invkur: array [0..2] of real;
begin
  whimres.descriptor:= 'D';
  setlength (whimres.value, 1);
  invkur[0]:= 0.0;
  invkur[1]:= 0.0;
  invkur[2]:= 0.0;
  if (sqs[0] <> 0) then begin invkur[0]:= (sqr(eigval[1]) * nat) / sqs[0]; end;
  if (sqs[1] <> 0) then begin invkur[1]:= (sqr(eigval[2]) * nat) / sqs[1]; end;
  if (sqs[2] <> 0) then begin invkur[2]:= (sqr(eigval[3]) * nat) / sqs[2]; end;
  whimres.value[0]:= invkur[0] + invkur[1] + invkur[2];
  whimDw:= whimres;
end;

function whimGkw (pcascores: TPCAScores; nat: integer; symtol: real): TWHIMresult;
var
  whimres: TWHIMresult;
  ns: real = 0;
  na: real = 0;
  i, k, m, j: integer;
  found: boolean = false;
begin
  whimres.descriptor:= 'G';
  setlength (whimres.value, 3);
  whimres.value[0]:= 0;
  whimres.value[1]:= 0;
  whimres.value[2]:= 0;
  for i:= 0 to 2 do begin
    na:= 0;
    ns:= 0;
    for j:= 0 to length(pcascores)-1 do begin
      found:= false;
      for k:= 0  to length(pcascores) - 1 do begin
      	if (j = k) then begin break; end;
        if (abs(pcascores[j][i] + pcascores[k][i]) < symtol) then begin
          found:= true;
        	ns:= ns + 1;
          break;
        end;
      end;
      if (found = false) then begin
        na:= na + 1;
      end;
    end;
	end;
  whimGkw:= whimres;
end;

function whimGw (pcascores: TPCAScores; nat: integer; symtol: real): TWHIMresult;
var
  whimres: TWHIMresult;
  eigvals: real;
begin
 	whimres.descriptor:= 'G';
  setlength (whimres.value, 1);
  whimres.value[0]:= 0;
end;

end.


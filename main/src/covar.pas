unit CoVar;

{$mode objfpc}{$H+}

interface

uses
Classes, SysUtils, Coord, Math, Typ, TypeDef;

const
  debye = 0.20819434;

type
  TWeights = array of TVector;
  TCovariance = array[1..3, 1..3] of arbFloat;
  TWCovariance = record
  	covariance: array of TCovariance;
  end;

function centralize (data: TVector): TVector;
function centralizeXYZ (data: TVVector): TVVector;
function centralizeXYZW (data: TVVector; w: TVector): TVVector;
function normalize (data: TVector; scale: real): TVector;
function standarize (data: TVector): TVector;
function getCovariance (coor: TVVector; w: TVector): TCovariance;



implementation

function scale (data: TVector; factor: real): TVector;
var
  value: real;
  std: TVector;
  i: integer = 0;
begin
  setLength (std, length(data));
	for value in data do begin
    std[i]:= value + factor;
    i:= i + 1;
  end;
  scale:= std;
end;

function centralize (data: TVector): TVector;
var
  value: real;
  avg: real;
  std: TVector;
  ssd: real;
  minval, maxval, sca: real;
  i: integer = 0;
begin
	setLength (std, length(data));
  avg:= mean(data);
  i:= 0;
  for value in data do begin
    std[i]:= (value - avg);
    i:= i + 1;
  end;
  centralize:= std;
end;

function centralizeXYZ (data: TVVector): TVVector;
var
  value: TVector;
  std: TVVector;
  i: integer = 0;
  avgx: real = 0;
  avgy: real = 0;
  avgz: real = 0;
  sumx: real = 0;
  sumy: real = 0;
  sumz: real = 0;
  sumw: real = 0;
begin
  setLength (std, length(data), 3);
  for value in data do begin
			sumx:= sumx + (value[0]);
			sumy:= sumy + (value[1]);
 			sumz:= sumz + (value[2]);
			sumw:= sumw + 1;
    i:= i + 1;
  end;
 	avgx:= sumx / sumw;
  avgy:= sumy / sumw;
 	avgz:= sumz / sumw;
  i:= 0;
	for value in data do begin
    std[i][0]:= value[0] - avgx;
    std[i][1]:= value[1] - avgy;
    std[i][2]:= value[2] - avgz;
    i:= i + 1;
  end;
  centralizeXYZ:= std;
end;

function centralizeXYZW (data: TVVector; w: TVector): TVVector;
var
  value: TVector;
  std: TVVector;
  i: integer = 0;
  avgx: real = 0;
  avgy: real = 0;
  avgz: real = 0;
  sumx: real = 0;
  sumy: real = 0;
  sumz: real = 0;
  sumw: real = 0;
begin
  setLength (std, length(data), 3);
  for value in data do begin
			sumx:= sumx + (value[0] * w[i]);
			sumy:= sumy + (value[1] * w[i]);
 			sumz:= sumz + (value[2] * w[i]);
			sumw:= sumw + w[i];
    i:= i + 1;
  end;
  if (sumw > 0) then begin
 		avgx:= sumx / sumw;
  	avgy:= sumy / sumw;
 		avgz:= sumz / sumw;
  end;
  i:= 0;
	for value in data do begin
    std[i][0]:= value[0] - avgx;
    std[i][1]:= value[1] - avgy;
    std[i][2]:= value[2] - avgz;
    i:= i + 1;
  end;
  centralizeXYZW:= std;
end;

function standarize (data: TVector): TVector;
var
  value: real;
  avg: real = 0;
  ssd: real = 0;
  std: TVector;
  i: integer = 0;
begin
  setLength (std, length(data));
  avg:= mean(data);
  ssd:= stddev(data);
  i:= 0;
  for value in data do begin
    std[i]:= ((value) - avg) / ssd;
    i:= i + 1;
  end;
  standarize:= std;
end;

function normalize (data: TVector; scale: real): TVector;
var
  value: real;
  min: real;
  max: real;
  nor: TVector;
  i: integer = 0;
begin
  min:= minValue (data);
  max:= maxValue (data);
  setLength (nor, length(data));
  if ((max-min) <> 0) then begin
   for value in data do begin
   	nor[i]:= scale + ((value - min) / (max - min));
   	i:= i + 1;
   end;
  	normalize:= nor;
  end else begin
  	normalize:= data;
  end;
end;

function getCovariance (coor: TVVector; w: TVector): TCovariance;
var
  value: TVector;
  weights: real = 0;
  x_avg: real = 0;
  y_avg: real = 0;
  z_avg: real = 0;
  ax: TVector;
  ay: TVector;
  az: TVector;
  xx_cov: real = 0;
  yy_cov: real = 0;
  zz_cov: real = 0;
  xy_cov: real = 0;
  xz_cov: real = 0;
  yz_cov: real = 0;
  cov: TCovariance;
  residuals: array of array[0..2] of real;
  i: integer = 0;
  k: integer = 0;
begin
  weights:= 0;
	i:= 0;
	for value in coor do begin
		xx_cov:= xx_cov + (value[0] * value[0] * w[i]);
		yy_cov:= yy_cov + (value[1] * value[1] * w[i]);
		zz_cov:= zz_cov + (value[2] * value[2] * w[i]);
		xy_cov:= xy_cov + (value[0] * value[1] * w[i]);
		xz_cov:= xz_cov + (value[0] * value[2] * w[i]);
		yz_cov:= yz_cov + (value[1] * value[2] * w[i]);
    //write (value[0], '  ', value[1],'  ',value[2]);
    //writeln;
    weights:= weights + w[i];
		i:= i + 1;
	end;
  if (weights <> 0) then begin
  	xx_cov:= xx_cov / weights;
		yy_cov:= yy_cov / weights;
		zz_cov:= zz_cov / weights;
		xy_cov:= xy_cov / weights;
		xz_cov:= xz_cov / weights;
		yz_cov:= yz_cov / weights;
  end else begin
  	xx_cov:= 0;
		yy_cov:= 0;
		zz_cov:= 0;
		xy_cov:= 0;
		xz_cov:= 0;
		yz_cov:= 0;
  end;
	cov[1, 1]:= xx_cov;
	cov[1, 2]:= xy_cov;
	cov[1, 3]:= xz_cov;
	cov[2, 1]:= xy_cov;
	cov[2, 2]:= yy_cov;
	cov[2, 3]:= yz_cov;
	cov[3, 1]:= xz_cov;
	cov[3, 2]:= yz_cov;
	cov[3, 3]:= zz_cov;
  //writeln (FloatToStr(cov[1, 1]), ' ', FloatToStr(cov[1, 2]), ' ', FloatToStr(cov[1, 3]));
 // writeln (FloatToStr(cov[2, 1]), ' ', FloatToStr(cov[2, 2]), ' ', FloatToStr(cov[2, 3]));
  //writeln (FloatToStr(cov[3, 1]), ' ', FloatToStr(cov[3, 2]), ' ', FloatToStr(cov[3, 3]));
  getCovariance:= cov;
end;

end.


unit Moments;

{$mode objfpc}{$H+}

interface

uses
Classes, SysUtils, CoVar, Typ, Math, Coord, TypeDef, ComProF;


type
	TMoments = record
   	descriptor: array [0..8] of string;     //0 shiftx 1 shifty 2 shiftz
   	value: array[0..8] of real;             //3 momxx, 4 momyy, 5 momzz, 6 momxy, 7 momyz
  end;
  TPolarMoments = record
  	descriptor: array [0..16] of string;  	//0 countp 1 countn 2 spol+ 3 spol- 4 dipolex 5 dipoley
    value: array [0..16] of real;           //6 dipolez 7 magxy 8 magxz 9 magyz 10 magxyz 11 dipxpo
  end;																			//12 dipypo 13 dipzpo 14 dipxno 15 dipyno 16 dipzno

function getMoment (coor: TCoordinates; w: TVector): TMoments;
function getPolarMoment (coor: TCoordinates; w: TVector): TPolarMoments;

implementation


function getMoment (coor: TCoordinates; w: TVector): TMoments;
var
  i: integer = 0;
  k: integer = 0;
  x: integer = 0;
  value: TVector;
  weights: real;
  poncoor: array of array [0..2] of real;
  geocenter: array [0..2] of real;
  polcenter: array [0..2] of real;
  shift: array [0..2] of real;
  momx: real = 0;
	momy: real = 0;
  momz: real = 0;
  momxx: real = 0;
	momyy: real = 0;
	momzz: real = 0;
	momxy: real = 0;
	momxz: real = 0;
	momyz: real = 0;
  moment: TMoments;
begin
  moment.descriptor[0]:= 'SHx';
	moment.descriptor[1]:= 'SHy';
	moment.descriptor[2]:= 'SHz';
  moment.descriptor[3]:= 'MOMxx';
	moment.descriptor[4]:= 'MOMyy';
	moment.descriptor[5]:= 'MOMzz';
	moment.descriptor[6]:= 'MOMxy';
	moment.descriptor[7]:= 'MOMxz';
  moment.descriptor[8]:= 'MOMyz';
  moment.value[0]:= 0.0;
	moment.value[1]:= 0.0;
	moment.value[2]:= 0.0;
	moment.value[3]:= 0.0;
	moment.value[4]:= 0.0;
	moment.value[5]:= 0.0;
	moment.value[6]:= 0.0;
	moment.value[7]:= 0.0;
	moment.value[8]:= 0.0;
  if ((((minvalue(w) >= 0) and (maxvalue(w) > 0)) or ((minvalue(w) < 0) and (maxvalue(w) <= 0))) and (minvalue(w) <> maxvalue(w))) then begin
    geocenter[0]:= 0;
		geocenter[1]:= 0;
		geocenter[2]:= 0;
		polcenter[0]:= 0;
		polcenter[1]:= 0;
		polcenter[2]:= 0;
		setlength (poncoor, length (coor.atom));
		for value in coor.atom do begin
			geocenter[0]:= geocenter[0] + value[0];
			geocenter[1]:= geocenter[1] + value[1];
			geocenter[2]:= geocenter[2] + value[2];
			polcenter[0]:= polcenter[0] + (value[0] * w[i]);
			polcenter[1]:= polcenter[1] + (value[1] * w[i]);
			polcenter[2]:= polcenter[2] + (value[2] * w[i]);
 			weights:= weights + w[i];
			i:= i + 1;
		end;
		for k:= 0 to 2 do begin
      geocenter[k]:= geocenter[k] / length (coor.atom);
  		polcenter[k]:= polcenter[k] / weights;
		end;
		shift[0]:= polcenter[0] - geocenter[0];
		shift[1]:= polcenter[1] - geocenter[1];
		shift[2]:= polcenter[2] - geocenter[2];
		i:= 0;
    momxx:= 0.0;
    momyy:= 0.0;
		momzz:= 0.0;
		momxy:= 0.0;
		momxz:= 0.0;
		momyz:= 0.0;
		for value in coor.atom do begin
      momxx:= momxx + (sqr(value[1]-polcenter[1]) + sqr(value[2]-polcenter[2])) * w[i];
      momyy:= momyy + (sqr(value[0]-polcenter[0]) + sqr(value[2]-polcenter[2])) * w[i];
			momzz:= momzz + (sqr(value[0]-polcenter[0]) + sqr(value[1]-polcenter[1])) * w[i];
			momxy:= momxy + ((value[0]-polcenter[0]) * (value[1]-polcenter[1])) * w[i];
			momxz:= momxz + ((value[0]-polcenter[0]) * (value[2]-polcenter[2])) * w[i];
			momyz:= momyz + ((value[1]-polcenter[1]) * (value[2]-polcenter[2])) * w[i];
			i:= i + 1;
		end;
		moment.descriptor[0]:= 'SHx';
		moment.descriptor[1]:= 'SHy';
		moment.descriptor[2]:= 'SHz';
		moment.descriptor[3]:= 'MOMxx';
		moment.descriptor[4]:= 'MOMyy';
		moment.descriptor[5]:= 'MOMzz';
		moment.descriptor[6]:= 'MOMxy';
		moment.descriptor[7]:= 'MOMxz';
    moment.descriptor[8]:= 'MOMyz';
		moment.value[0]:= roundd(shift[0], decimal);
		moment.value[1]:= roundd(shift[1], decimal);
		moment.value[2]:= roundd(shift[2], decimal);
		moment.value[3]:= roundd(momxx, decimal);
		moment.value[4]:= roundd(momyy, decimal);
		moment.value[5]:= roundd(momzz, decimal);
		moment.value[6]:= roundd(momxy, decimal);
		moment.value[7]:= roundd(momxz, decimal);
		moment.value[8]:= roundd(momyz, decimal);
	end;
	getMoment:= moment;
end;


function getPolarMoment (coor: TCoordinates; w: TVector): TPolarMoments;
var
  i: integer = 0;
  k: integer = 0;
  z: integer = 0;
  t: integer = 0;
  x: integer = 0;
  countp: integer = 0;
  countn: integer = 0;
  value: TVector;
  magnitude: TpolarMoments;
  spol: array [0..1] of real;
  poncoor: array of array [0..2] of real;
  momcoor: array of array [0..2] of real;
  geocenter: array [0..2] of real;
  polcenter: array [0..1] of array [0..2] of real;
  dipole: array [0..2] of real;
	magxy: real = 0;
	magxz: real = 0;
	magyz: real = 0;
	magxyz: real = 0;
	dipxpo: real = 0;
	dipypo: real = 0;
	dipzpo: real = 0;
	dipxno: real = 0;
	dipyno: real = 0;
	dipzno: real = 0;
  d: real = 0;
  h: real = 0;
begin
  magnitude.descriptor[0]:= 'CTp';
	magnitude.descriptor[1]:= 'CTn';
	magnitude.descriptor[2]:= 'QTp';
	magnitude.descriptor[3]:= 'QTn';
	magnitude.descriptor[4]:= 'Dx';
	magnitude.descriptor[5]:= 'Dy';
	magnitude.descriptor[6]:= 'Dz';
	magnitude.descriptor[7]:= 'Mxy';
	magnitude.descriptor[8]:= 'Mxz';
	magnitude.descriptor[9]:= 'Myz';
	magnitude.descriptor[10]:= 'Mxyz';
	magnitude.descriptor[11]:= 'SHpxo';
	magnitude.descriptor[12]:= 'SHpyo';
	magnitude.descriptor[13]:= 'SHpzo';
	magnitude.descriptor[14]:= 'SHnxo';
	magnitude.descriptor[15]:= 'SHnyo';
	magnitude.descriptor[16]:= 'SHnzo';
  magnitude.value[0]:= 0;
	magnitude.value[1]:= 0;
	magnitude.value[2]:= 0;
	magnitude.value[3]:= 0;
	magnitude.value[4]:= 0;
	magnitude.value[5]:= 0;
	magnitude.value[6]:= 0;
	magnitude.value[7]:= 0;
	magnitude.value[8]:= 0;
	magnitude.value[9]:= 0;
	magnitude.value[10]:= 0;
	magnitude.value[11]:= 0;
	magnitude.value[12]:= 0;
	magnitude.value[13]:= 0;
	magnitude.value[14]:= 0;
	magnitude.value[15]:= 0;
	magnitude.value[16]:= 0;
  if ((minvalue(w) < 0) and (maxvalue(w) > 0)) then begin
  	setlength (poncoor, length (coor.atom));
		polcenter[0][0]:= 0.0;
		polcenter[0][1]:= 0.0;
		polcenter[0][2]:= 0.0;
		polcenter[1][0]:= 0.0;
		polcenter[1][1]:= 0.0;
		polcenter[1][2]:= 0.0;
		for value in coor.atom do begin
			geocenter[0]:= geocenter[0] + value[0];
			geocenter[1]:= geocenter[1] + value[1];
			geocenter[2]:= geocenter[2] + value[2];
			if (w[i] > 0) then begin
				polcenter[0][0]:= polcenter[0][0] + (value[0] * w[i]);
				polcenter[0][1]:= polcenter[0][1] + (value[1] * w[i]);
				polcenter[0][2]:= polcenter[0][2] + (value[2] * w[i]);
				spol[0]:= spol[0] + w[i];
				countp:= countp + 1;
			end;
			if (w[i] < 0) then begin
				polcenter[1][0]:= polcenter[1][0] - (value[0] * w[i]);
				polcenter[1][1]:= polcenter[1][1] - (value[1] * w[i]);
				polcenter[1][2]:= polcenter[1][2] - (value[2] * w[i]);
				spol[1]:= spol[1] - w[i];
				countn:= countn + 1;
			end;
			i:= i + 1;
		end;
		geocenter[0]:= geocenter[0] / length (coor.atom);
		geocenter[1]:= geocenter[1] / length (coor.atom);
		geocenter[2]:= geocenter[2] / length (coor.atom);
		for t:= 0 to 1 do begin
				for k:= 0 to 2 do begin
					polcenter[t][k]:= polcenter[t][k] / spol[t];
				end;
		end;
		dipole[0]:= polcenter[0][0] - polcenter[1][0];
		dipole[1]:= polcenter[0][1] - polcenter[1][1];
		dipole[2]:= polcenter[0][2] - polcenter[1][2];
		magxy:= sqrt(sqr(dipole[0]) + sqr(dipole[1]));
		magxz:= sqrt(sqr(dipole[0]) + sqr(dipole[2]));
		magyz:= sqrt(sqr(dipole[1]) + sqr(dipole[2]));
		magxyz:= sqrt(sqr(dipole[0]) + sqr(dipole[1]) + sqr(dipole[2]));
		dipxpo:= polcenter[0][0]-geocenter[0];
		dipypo:= polcenter[0][1]-geocenter[1];
		dipzpo:= polcenter[0][2]-geocenter[2];
		dipxno:= polcenter[1][0]-geocenter[0];
		dipyno:= polcenter[1][1]-geocenter[1];
		dipzno:= polcenter[1][2]-geocenter[2];
		magnitude.value[0]:= roundd(countp, decimal);
		magnitude.value[1]:= roundd(countn, decimal);
		magnitude.value[2]:= roundd(spol[0], decimal);
		magnitude.value[3]:= roundd(spol[1], decimal);
		magnitude.value[4]:= roundd(dipole[0], decimal);
		magnitude.value[5]:= roundd(dipole[1], decimal);
		magnitude.value[6]:= roundd(dipole[2], decimal);
		magnitude.value[7]:= roundd(magxy, decimal);
		magnitude.value[8]:= roundd(magxz, decimal);
		magnitude.value[9]:= roundd(magyz, decimal);
		magnitude.value[10]:= roundd(magxyz, decimal);
		magnitude.value[11]:= roundd(dipxpo, decimal);
		magnitude.value[12]:= roundd(dipypo, decimal);
		magnitude.value[13]:= roundd(dipzpo, decimal);
		magnitude.value[14]:= roundd(dipxno, decimal);
		magnitude.value[15]:= roundd(dipyno, decimal);
		magnitude.value[16]:= roundd(dipzno, decimal);
  end;
  getPolarMoment:= magnitude;
end;


end.


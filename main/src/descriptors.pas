unit Descriptors;

{$mode objfpc}{$H+}

interface

uses
  SysUtils, Math, Typ, Eig, Omv, Mol2, STDMol, OFileTextArray, CoVar, Coord, TypeDef,
  Weights, WHIM, Moments, ComProF;

const
  WHIMdescriptors = 8;

function csvDescriptor (fmolname: string; fparname: string): TDataCSV;

implementation

function csvDescriptor (fmolname: string; fparname: string): TDataCSV;
var
  whimres: array of array of TWHIMresult;     				//Matriz con los resultados WHIM
  moment: array of TMoments;                     			//Matriz de momentos por cada esquema de pesada
  polmag: array of TPolarMoments;                     //Matriz de momentos por cada esquema de pesada;                     //Matriz de momentos por cada esquema de pesada
  comCoordinates: TCoordinates;    										//Registro con coordenadas atomicas compiladas
  comCenterCoorM: array of TVVector;                	//Matriz centralizada de coordenadas
  comCenterCoorMT: array of array of array of real; 	//Matriz centralizada de coordenadas transpuesta
  comCenterCoorMW: array of TVVector;               	//Matriz centralizada de coordenadas pesadas
  comCenterCoorMWT: array of array of array of real;  //Matriz centralizada de coordenadas pesadas transpuesta
  comWeights: TWcompiled;                         		//Registro con pesos por cada atomo compilados
  comWeightsNor: array of TVector;                    //Registro con pesos por cada atomo compilados
  molecule: TTripos;                              		//Registro de molecula con estructura TRIPOS
  covariance: array of TCovariance;               		//Matriz de covarianzas por cada esquema de pesada
 	xchgeval: real;                                     //Variable de intercambio al reordenar los autovalores
  xchgevec: array[1..3] of ArbFloat;                  //Matriz de intercambio al reordenar los autovectors
  eigval: array of TWHIMEigenvalues;                  //Matriz de N=NPesos con su autovalor
  eigvec: array of TWHIMEigenvectors;                 //Matriz de N=NPesos con su autovector
  eigvecM: TWHIMEigenvectors;                         //Variable de trabajo
  eigvecMT: TWHIMEigenvectors;                        //Variable de trabajo
  eigenvalues: TEigenvalues;                          //Variable de trabajo
  eigenvectors: TEigenvectors;                        //Variable de trabajo
  pcaScores: array of TPCAScores;                     //Matriz de 3xNatomos con los scores PCA
  pcaScoresT: array of TPCAScores;                    //Su transpuesta
  pcaSqscores: array [0..2] of real;                  //Matriz con la suma de los scores^4
  pcaSqscore: real = 0;                               //Sumatoria de scores^4
  lines: array of string;												  		//Matriz que contiene el archivo mol2 por linea
  term: integer = 0;                                  //Retorno de error para el calculo de autovectores/valores
  x, r, i, k, z: integer;                             //Variables de uso general
  wlist: TWeightList;                                 //Matriz de pesos: (archivo parametros) + (molumna q mol2)
  csv: TDataCSV;                                      //Estructura CSV a devolver
begin

  loadWeights (fparname);
  wlist:= getWeights();
  lines:= readFileText (fmolname);
  molecule:= mol2Process (lines);
  setLength (wlist, length(wlist) + 1);
  setLength (covariance, length(wlist));
  setLength (comWeights, length(wlist), molecule.nat);
  setLength (comCoordinates.atom, molecule.nat, 3);
	setLength (comCoordinates.weight, molecule.nat);
  setLength (comCoordinates.atype, molecule.nat);
  setLength (comCoordinates.charge, molecule.nat);
  wlist[length(wlist)-1]:= 'q';
  for i:= 0 to (molecule.nat - 1) do begin
    comCoordinates.atom[i,0]:= molecule.atoms[i].x;
    comCoordinates.atom[i,1]:= molecule.atoms[i].y;
    comCoordinates.atom[i,2]:= molecule.atoms[i].z;
    comCoordinates.atype[i]:= molecule.atoms[i].atom_type;
    comCoordinates.charge[i]:= molecule.atoms[i].charge;
  end;
  for r:= 0 to length(wlist) - 1 do begin
    for i:= 0 to (molecule.nat - 1) do begin
      	if (r = (length(wlist) - 1)) then begin
        	comWeights[r][i]:= molecule.atoms[i].charge;
        end else begin
          comWeights[r][i]:= getWeightvalue(comCoordinates.atype[i], r);

        end;
    end;
  end;
  setLength (comCenterCoorM, length (covariance), molecule.nat, 3);
 	setLength (comCenterCoorMT, length (covariance), 3, molecule.nat);
  setLength (comCenterCoorMW, length (covariance), molecule.nat, 3);
 	setLength (comCenterCoorMWT, length (covariance), 3, molecule.nat);
  setLength (eigval, length (covariance));
  setLength (eigvec, length (covariance));
  setLength (whimres, length (covariance));
  setlength (pcaScores, length (covariance), 3, molecule.nat);
  setlength (pcaScoresT, length (covariance), molecule.nat, 3);
  setLength (comWeightsNor, length (covariance), molecule.nat);
  setLength (moment, length (covariance));
  setLength (polmag, length (covariance));
  for r:= 0 to length (covariance) - 1 do begin
   	if (minvalue(comWeights[r]) < 0) then begin
      comWeightsNor[r]:= normalize (comWeights[r], 1);
    end else begin
      comWeightsNor[r]:= comWeights[r];
    end;
  	for i:= 0 to (molecule.nat - 1) do begin
    	comCenterCoorM[r][i][0]:= molecule.atoms[i].x;
    	comCenterCoorM[r][i][1]:= molecule.atoms[i].y;
    	comCenterCoorM[r][i][2]:= molecule.atoms[i].z;
    end;
    comCenterCoorM[r]:= centralizeXYZ (comCenterCoorM[r]);
    comCenterCoorMW[r]:= centralizeXYZW (comCenterCoorM[r], comWeightsNor[r]);
    for k:= 0 to 2 do	begin
  		for z:= 0 to molecule.nat - 1 do begin
        comCenterCoorMT[r][k][z]:= comCenterCoorM[r][z][k];
        comCenterCoorMWT[r][k][z]:= comCenterCoorMW[r][z][k];
      end;
  	end;
  end;
  for r:= 0 to length (covariance) - 1 do begin
    covariance[r]:= getCovariance(comCenterCoorM[r], comWeightsNor[r]);
    eigge3(covariance[r][1,1], 3, 3, eigenvalues[1], eigenvectors[1,1], 3, term);
    eigval[r][1]:= eigenvalues[1].re;
    eigval[r][2]:= eigenvalues[2].re;
    eigval[r][3]:= eigenvalues[3].re;
    eigvec[r][1][1]:= (eigenvectors[1,1].re * 1);
    eigvec[r][1][2]:= (eigenvectors[2,1].re * 1);
    eigvec[r][1][3]:= (eigenvectors[3,1].re * 1);
    eigvec[r][2][1]:= (eigenvectors[1,2].re * 1);
    eigvec[r][2][2]:= (eigenvectors[2,2].re * 1);
    eigvec[r][2][3]:= (eigenvectors[3,2].re * 1);
    eigvec[r][3][1]:= (eigenvectors[1,3].re * 1);
    eigvec[r][3][2]:= (eigenvectors[2,3].re * 1);
    eigvec[r][3][3]:= (eigenvectors[3,3].re * 1);
  	for z:= 1 to 3 do	begin
      for i:= 1 to 3 do begin
        if ((eigval[r][i]) < (eigval[r][z])) then begin
          xchgeval:= eigval[r][z];
          xchgevec:= eigvec[r][z];
          eigval[r][z]:= eigval[r][i];
          eigval[r][i]:= xchgeval;
          eigvec[r][z]:= eigvec[r][i];
         	eigvec[r][i]:= xchgevec;
        end;
      end;
    end;
    for z:= 1 to 3 do	begin
      for i:= 1 to 3 do begin
        eigvecM[z][i]:= eigvec[r][z][i];
      end;
    end;
    omvtrm(eigvecM[1,1], 3, 3, 3, eigvecMT[1,1], 3);
    for z:= 0 to 2 do begin
    	for i:= 0 to molecule.nat - 1 do begin
        	for k:=0 to 2 do begin
            pcaScores [r][z,i]:= pcaScores [r][z,i] + (eigvecM[z+1, k+1] * comCenterCoorMT[r][k, i]);
          end;
       end;

  	end;
    for k:= 0 to 2 do	begin
  		for z:= 0 to molecule.nat - 1 do begin
        pcaScoresT[r][z][k]:= pcaScores[r][k][z];
      end;
  	end;
   	for z:=0 to 2 do begin
      pcaSqscore:= 0;
      for k:=0 to molecule.nat do begin
        pcaSqscore:= pcaSqscore + sqr(sqr((pcaScores[r][z][k])));
      end;
       pcaSqscores[z]:= pcaSqscore;
   	end;
   	for x:= 0 to WHIMdescriptors - 1 do begin
       setLength (whimres[r], WHIMdescriptors);
       case x of
       		0: begin whimres[r][x]:= whimLkw (eigval[r]); end;
       		1: begin whimres[r][x]:= whimPkw (eigval[r]); end;
        	2: begin whimres[r][x]:= whimEkw (eigval[r], pcaSqscores, molecule.nat); end;
       		3: begin whimres[r][x]:= whimTw (eigval[r]); end;
        	4: begin whimres[r][x]:= whimAw (eigval[r]); end;
       		5: begin whimres[r][x]:= whimVw (eigval[r]); end;
        	6: begin whimres[r][x]:= whimKw (eigval[r]); end;
       		7: begin whimres[r][x]:= whimDw (eigval[r], pcaSqscores, molecule.nat); end;
        	8: begin whimres[r][x]:= whimGkw (pcaScoresT[r], molecule.nat, 0.01); end;
        	9: begin whimres[r][x]:= whimGw (pcaScores[r], molecule.nat, 0.01); end;
       end;
    end;
    moment[r]:= getMoment (comCoordinates, comWeights[r]);
    polmag[r]:= getPolarMoment (comCoordinates, comWeights[r]);
 	end;
  csv.fields:= 0;
  for k:=0 to length (whimres) - 1 do begin
    for x:= 0 to length(whimres[k]) - 1 do begin
    	for i:= 0 to length(whimres[k][x].value) - 1 do	begin
				setlength (csv.field, length(csv.field)+1);
        setlength (csv.value, length(csv.value)+1);
        if (length(whimres[k][x].value) > 1) then begin
       		csv.field[csv.fields]:= whimres[k][x].descriptor + inttostr(i+1) + wlist[k];
        end else begin
          csv.field[csv.fields]:= whimres[k][x].descriptor + wlist[k];
        end;
        csv.value[csv.fields]:= floattostr(roundd(whimres[k][x].value[i], decimal));
        csv.fields:= csv.fields + 1;
       	end;
      end;
  end;
 	for k:=0 to length (moment) - 1 do begin
    for x:= 0 to length(moment[k].descriptor) - 1 do begin
    	setlength (csv.field, length(csv.field)+1);
      setlength (csv.value, length(csv.value)+1);
      csv.field[csv.fields]:= moment[k].descriptor[x] + wlist[k];
      csv.value[csv.fields]:= floattostr(roundd(moment[k].value[x], decimal));
      csv.fields:= csv.fields + 1;
    end;
  end;
  for k:=0 to length (polmag) - 1 do begin
    for x:= 0 to length(polmag[k].descriptor) - 1 do begin
    	setlength (csv.field, length(csv.field)+1);
      setlength (csv.value, length(csv.value)+1);
      csv.field[csv.fields]:= polmag[k].descriptor[x] + wlist[k];
      csv.value[csv.fields]:= floattostr(roundd(polmag[k].value[x], decimal));
      csv.fields:= csv.fields + 1;
    end;
  end;
  csvDescriptor:= csv;
end;

end.


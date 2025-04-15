program QTDD;

uses
  SysUtils, Math, Typ, Eig, Omv, Mol2, STDMol, OFileTextArray, CoVar, Coord, TypeDef,
  Weights, WHIM, Moments, Descriptors, Dos, LazUtils, comprof;

{$I+}

var
  csv: TDataCSV;
  name: string;
  parf: string = '';
  path: string = '';
  wdir: string = './';
  outf: string = '';
  outfh: text;
  mask: string = '*.mol2';
  drec: SearchRec;
  csvsep: string = '';
  i, k: integer;
  biof: string = '';
  csvline: string;
  csvlnid: longint;
  lnerror: longint;
  hour, min, sec, hsec : word;
  year, mnth, day, wday : word;

begin
  writeln;
  writeln (' ATOQ Amber To QSAR');
  writeln (' 2020 Facultad de Farmacia y Bioquimica, Universidad de Buenos Aires.');
  writeln (' 2020 Gabriel Jasinski, Lucas Fabian, Florencia Martini, Alberina Moglioni.');
  writeln;

  k:= 0;
  i:= 0;
  outf:= '';
  biof:='';
  wdir:= '.';
  parf:= './atoq.par.csv';
  write ('+ Command line options: ');
  while k <= paramcount + 1 do begin
    k:= k + 1;
    i:= i + 1;
    write (paramstr(i), ' ');
   	if ((paramstr(k) = '-bio') and (k < paramcount)) then begin
     biof:= (paramstr(k + 1));
     k:= k + 1;
    end;
    if ((paramstr(k) = '-in') and (k < paramcount)) then begin
     wdir:= (paramstr(k + 1));
     k:= k + 1;
    end;
    if ((paramstr(k) = '-out') and (k < paramcount)) then begin
     outf:= (paramstr(k + 1));
     k:= k + 1;
    end;
    if ((paramstr(k) = '-par') and (k < paramcount)) then begin
     parf:= (paramstr(k + 1));
     k:= k + 1;
    end;
  end;

  writeln;
  writeln;

  if (length(outf) = 0) then begin
   gettime(hour, min, sec, hsec);
   getdate(year, mnth, day, wday);
   outf:= wdir + '/' + inttostr(year) + inttostr(mnth) + inttostr(day) + inttostr(hour) + inttostr(min) + '.desc.csv';
  end;

  if (directoryexists(wdir) = false) then begin
    writeln ('## ERROR: Working directory ', wdir, ' does not exist.');
    writeln;
    writeln ('Press ENTER to exit');
    readln;
    exit;
  end;

  if (fileexists(parf) = false) then begin
   	writeln ('## ERROR: Parameter file ', parf, ' does not exist.');
    writeln;
    writeln ('Press ENTER to exit');
    readln;
    exit;
  end;

 	try
 		assign (outfh, outf);
   	rewrite (outfh);
 	except
  	writeln ('## ERROR: The path of file ', outf, '  is incorrect.');
    writeln;
		writeln ('Press ENTER to exit');
		readln;
    close(outfh);
		exit;
 	end;

  writeln ('+ Notice: Normalization with scaling of the polar weights will be used to calculate WHIM');
 	writeln ('+ Notice: The molecular moments will be calculated with only the non-polar weights');
	writeln ('+ Notice: The molecular dipole magnitude will be calculated with only the polar weights');
  writeln;
  writeln ('- Input directory: ', wdir);
  writeln ('- Parameter file:  ', parf);
  writeln ('- Output file:     ', outf);
  writeln;

 	csvlnid:= 0;
  lnerror:= 0;
  findfirst(wdir + '/' + mask, archive, drec);
	while (DOSerror = 0) do begin
  	path:= wdir + '/' + drec.name;
  	if (fileexists (path)) then begin
      name:= split(drec.name, '.')[0];
      write ('- Molecule: ', name, ' ');
      try
    		csv:= csvDescriptor (path, parf);
        csvsep:= ',';
     		if (csvlnid = 0) then begin
      		csvline:= 'File,';
     			for k:= 0 to csv.fields - 1 do begin
      			if (k = csv.fields - 1) then begin  csvsep:= ''; end;
          	csvline:= csvline + csv.field[k] + csvsep;
     			end;
        	writeln (outfh, csvline);
     		end;
      	csvline:= '"'+ name + '",';
      	csvsep:= ',';
  	 		for k:= 0 to csv.fields - 1 do begin
      		if (k = csv.fields - 1) then begin  csvsep:= ''; end;
        	csvline:= csvline + csv.value[k] + csvsep;
     		end;
     		writeln (outfh, csvline);
     		csvlnid:= csvlnid + 1;
        write ('OK');
      except
        	on E: exception do begin write (E.Message); lnerror:= lnerror + 1; end;
      end;
  	end;
    writeln;
    findnext(drec);
  end;
  writeln;
  findclose(drec);
  close (outfh);
  writeln ('+ ', (csvlnid), ' Processed files.');
  writeln ('+ ', (lnerror), ' Error files.');
  writeln;
  writeln ('Press ENTER to exit');
  readln;
end.




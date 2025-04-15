unit OFileTextArray;

{$mode objfpc}{$H+}

interface
uses
Classes, SysUtils;

type
	TTextArray = array of AnsiString;

function readFileText (fname: string): TTextArray;

implementation

function readFileText (fname: string): TTextArray;
var
	line_index: integer = 0;
  file_id: textfile;
  lines: TTextArray;
begin
  {$I-}
  assign (file_id, fname);
  {$I+}
	reset (file_id);
	while not eof(file_id) do
	begin
    line_index:= line_index + 1;
    setLength(lines, line_index + 1);
		readLn (file_id, lines[line_index]);
	end;
  close (file_id);
  readFileText:= lines;
end;

end.


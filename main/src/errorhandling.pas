unit ErrorHandling;

{$mode objfpc}{$H+}

interface

uses
Classes, SysUtils;

type
  TError = record
    code: string;
    name: string;
  end;

var
	error: TError;

procedure CatchUnhandledException(Obj: TObject; Addr: Pointer; FrameCount: Longint; Frames: PPointer);

implementation

end.


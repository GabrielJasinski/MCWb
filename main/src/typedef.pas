unit TypeDef;

{$mode objfpc}{$H+}

interface

uses
Classes, SysUtils, Math, Typ;

const
  decimal = 5;

type
  TVector = array of real;
  TSVector = array of string;
  TVVector = array of array of real;
  TEigenvalues = array[1..3] of complex;
  TEigenvectors = array[1..3, 1..3] of complex;
  TDataCSV = record
    fields: integer;
    field: array of string;
    value: array of string;
  end;



implementation



end.


                            (* SNOWBALL.PAS *)

(* This is a program to compute estimates according to the methods
 in Frank & Snijders (1994).
 Required format of input file:
 Vertices are indicated by numbers in the range 1 ... vmax
 (vmax depending on implementation).
 Each line must have the vertex number and then the numbers of other vertices
 to which an arrow from this vertex points.
 Blanks after the last number on a line lead to input errors.
 If indegrees are also available, these can be given in a separate
 input file (format: each line has first a vertex number and then the
 corresponding indegree), and then the Horvitz-Thompson estimators also will be
 computed. *)


{$R+}
{$N-}
{$M 40000, 0, 655360}  (*  extra stack size *)

uses Crt,Dos;

const vmax = 1000; outmax = 40; nmax = 400;
(* vmax is maximum number of respondents + nominees;
   nmax is maximum number of respondents, i.e., maximum initial sample size;
   outmax is maximum out-degree. *)

type vec = array[1..vmax] of Boolean;
     veci = array[1..vmax] of word;
     vecr = array[1..vmax] of real;
     vecrn = array[1..nmax] of real;

var n,m: integer;
   s0, s1: veci;
   sample0, sample1: vec;
   zeros,a, out, inb: veci;
   outball: array[1..nmax,1..outmax] of integer;
   indegrees: boolean;
   fn1,fn2: string;
   fn: string;
   outf: text;

procedure initialise;
var i,j: integer;
  aa: char;
begin
  clrscr;
  writeln('                          SNOWBALL ');
  writeln('                     <<<<          >>>>   ');
  writeln('                          --------       ');
  writeln;
  writeln('First read the file <README.TXT> to know the required input file.');
  writeln('If you have not read the README file, then it is best to stop now');
  writeln('and read the file first.');
  writeln;
  write('Do you want to stop? (Y/N): ');
  aa := readkey; writeln(aa);
  if (aa in ['y','Y']) then Halt;
  writeln;
  write('Give name for outputfile (default = SNOWBALL.OUT) : ');
  readln(fn); if (fn = '') then fn := 'SNOWBALL.OUT';
  assign(outf,fn);
  rewrite(outf);
  writeln(outf,'Estimation of population size using the methods in');
  writeln(outf,'O. Frank & T.A.B. Snijders,');
  writeln(outf,'Estimating the size of hidden populations using snowball sampling,');
  writeln(outf,'Journal of Official Statistics, 10 (1994), 53 - 67.');
  writeln(outf,'Symbols used are as in this paper.');
  writeln(outf);
  writeln(outf,'Current implementation: ');
  writeln(outf,'Vertex numbers of respondents and nominees must range between 1 and ',vmax,'.');
  writeln(outf,'Initial sample size maximum ',nmax,'.');
  writeln(outf,'Out-degree maximum ',outmax,'.');
  writeln(outf,'Larger implementation limits can be obtained by changing');
  writeln(outf,
    'the corresponding constants in the Turbo Pascal source code and recompiling.');
  writeln(outf);

  for i:=1 to vmax do begin
    s0[i] := 0; s1[i] := 0;
    a[i] := 0; out[i] := 0; inb[i] := 1;
    zeros[i] := 0;
    sample0[i] := false; sample1[i] := false;
   end;
for i:=1 to nmax do for j:=1 to outmax do outball[i,j] := 0;

end; (* initialise *)

procedure Waitkey;
var a: char;
begin
  writeln; writeln('              (press a key ....)');
  a := readkey;
end;

procedure readfile;
var i,j,kk,dum: integer;
    inf1,inf2: Text;
begin
write('Give filename for input : ');
repeat readln(fn1) until (fn1 > '');
assign(inf1,fn1);
writeln('Each line must have the vertex number');
writeln(' and then the numbers of other vertices');
writeln(' to which an arrow from this vertex points.');
writeln('Blanks after the last number on a line lead to input errors.');

reset(inf1);
writeln;
writeln('Reading input file ',fn1,':');
writeln;
n := 0;
while not eof(inf1) do begin
  Read(inf1,i);
  if (i >= 1) then begin
   Inc(n);
   s0[n] := i;
   writeln('Reading vertex number ',n ,' coded as ',i,'. Arrows point to: ');
   j := 0;
   if (n > nmax) then begin
     writeln('Initial sample size n in this version may not be higher than ',n,'.');
     writeln('Adapt constant nmax in the program.');
     waitkey; close(inf1); Halt;
    end;
   if (i > nmax) then begin
     writeln('Initial vertex number in this version may not be higher than ',nmax,'.');
     writeln('Adapt dataset or constant nmax in the program.');
     waitkey; close(inf1); Halt;
    end;
   while not Eoln(inf1) do begin
     Inc(j);
     if (j > outmax) then begin
       writeln('Out-degree in this version may not be higher than ',outmax,'.');
       writeln('You reported too many arcs from vertex ',i,'.');
       writeln('Adapt dataset or constant outmax in the program.');
       waitkey; close(inf1); Halt;
      end;
     read(inf1,dum);
     write(dum,' ');
     if (dum > vmax) then begin
       writeln('Mentioned vertex number in this version may not be higher than ',outmax,'.');
       writeln('You reported an arc from vertex ',i,' to vertex ',dum,'.');
       writeln('Adapt dataset or constant outmax in the program.');
       waitkey; close(inf1); Halt;
      end;
     outball[i,j] := dum;
    end;
    a[n] := j; out[i] := j;
    writeln;
    sample0[i] := true;
   end; (* if i >= 1 *)
 end; (* while not eof *)
close(inf1);

(* compute first wave snowball sample *)
m := 0;
for i:=1 to n do begin
  for kk := out[i] downto 1 do begin
    j := outball[i,kk];
    if (not (sample1[j] or sample0[j])) then begin
      sample1[j] := true;
      Inc(m);
      s1[m] := j;
     end; (* if *)
   end; (* kk *)
 end; (* i *)
writeln;
writeln('Give a filename for in-degrees, if these are available;');
writeln('each line must have the vertex number and its indegree.');
writeln(' (Give a blank return if the in-degree file does not exist.)  ');
write('Filename : ');
readln(fn2);
indegrees := (fn2 <> '');
if indegrees then begin
  assign(inf2,fn2);
  reset(inf2);
  while not eof(inf2) do readln(inf2,i,inb[i]);
 end;
writeln(outf);
write(outf,'Data read from file ',fn1);
if (indegrees) then write(outf,' and ',fn2);
writeln(outf,'.');

writeln(outf,'Vertices in initial sample:');
write(outf,'number    code   outdegree ');
if indegrees then writeln(outf,'  indegree  ') else writeln(outf);
for i:=1 to vmax do if (s0[i] > 0) then if sample0[s0[i]] then begin
  write(outf,i:5,' ',s0[i]:8,' ',a[i]:7);
  if indegrees then writeln(outf,inb[s0[i]]:11) else writeln(outf);
 end; (*if sample0 *)

end; (* readfile *)


function di(a,b : real): real;
begin
  if (b > 1e-20) then di := a/b else di := 0.0;
end;

function d(a, b: real):real;
begin
  if (b >= 1e-20) then d := a/b else d := a;
end;

function power(a: real; b: byte): real;
var x: real; i: byte;
begin
  x := 1;
  for i := b downto 1 do x := x*a;
  power := x;
end;

function max(a,b: real):real;
begin
  if (a < b) then max := b else max := a;
end;

function root(a:real):real;
begin
  if (a > 1e-20) then root := sqrt(a) else root := 0.0;
end;

function answeryes: boolean;
var a: char;
begin
repeat a := readkey until (a in ['y','Y','n','N']);
answeryes := (a in ['y','Y']);
if (a in ['y','Y']) then write('Y ') else write('N ');
end;

function v3_ex(nn,mm: integer; tt: real; v00: real; note: boolean): real;
(* Exact solution of the definition formula for v_3;
   see top of p. 57.
   Newton-Raphson algorithm is used for solving this equation. *)
var n,m: integer;
    t,vv,rrr,vmin: real;
    conver: boolean;
(* first two procedures *)

function fv5: real;
var tmm,ff: real;
begin
 if (t - 0.1 < m) then tmm := 1 else tmm := t-m;
 ff := n + m*(2*n+t-2)/(2*tmm) - sqr(m)/(2*(2*n+m-2));
 if (ff < vmin) then fv5 := vmin else fv5 := ff;
end; (* fv5 *)

procedure find(var vv: real; var conv: boolean);
var vv0,rrr,vvold,f,fprime: real;
begin
  vv0 := vv;
  conv := true;
  repeat
   vvold := vv0;
   f := n*ln(1.0 - t/(n*(vv0-1.0))) - ln(1.0 - m/(vv0-n));
   fprime := n*t/((vv0-1)*(n*(vv0-1)-t)) - m/((vv0-n)*(vv0-n-m));
   if (abs(fprime) < 1e-12) then begin
     vv0 := fv5;
     vvold := vv0;
     if note then write('x');
     conv := false;
    end else
     vv0 := vv0 - f/fprime; (* Newton-Raphson step *)
   vv0 := max(vv0,vmin);
   if note then write('.');
  until (abs(vv0 - vvold) < 0.1);
  vv := vv0;
end; (* find *)

begin (* v3_ex *)
  n := nn; m := mm; t := tt;
(* First determine initial value. *)
  vmin := n + m + 1;
  if (t-0.1 < m) then begin
    v3_ex := vmin; exit; end;
  if (v00 < vmin) then vv := fv5 else vv := v00;
  vv := max(vv,vmin);
(* Now iterate to find the solution. *)
  find(vv,conver);
  if (not conver) (* no convergence *) then begin
  (* try again with a very small starting value *)
    vv := vmin;
    find(vv,conver);
   end;
  v3_ex := vv;
end; (* v3_ex *)

procedure jack(a1: vecrn; var s: real);
(* mean and standard error *)
var i: integer;
    m0,s0: real;
begin
  m0 := 0; s0 := 0;
  for i := 1 to n do begin
    m0 := m0 + a1[i];
    s0 := s0 + sqr(a1[i]);
   end;
  m0 := m0/n;
  s := root((n/2-1.0)*(s0/n - sqr(m0)));
end;

procedure calcvhats1(n1,m1: integer; k,s,r: real;
          var v1,v30,v31,v32,v33,v4,v5,se1,se30,se31,se32,se33: real);
var t,tmm,rrr,vmin: real;
procedure ma(var vv:real);
begin
  if (vv < vmin) then vv := vmin;
end;
begin
  vmin := n1+m1;
  t := r+s;
  if (t-0.1 < m1) then tmm := 1 else tmm := t-m1;
  v1 := n1 + (n1-1.0)*s/max(r,1);
  v30 := n1 + (n1-1.0)*m1/tmm;
  rrr := 8.0*(n1-1.0)*(t-m1)/(sqr(2.0*n1+m1-2.0));
  v31 := n1 + (m1*(2.0*n1+m1-2.0)*(1.0+sqrt(1+rrr)))/(4*tmm);
  v32 := n1 + m1*(2.0*n1+t-2)/(2*tmm) - sqr(1.0*m1)/(2.0*(2.0*n1+m1-2));
  v33 := n1 + m1*(2.0*n1+t-2)/(2*tmm);
  v4 := n1 + 1.0*n1*m1/max(k,1);
  v5 := n1 + (n1-1.0)*m1/max(k,1);
  ma(v1); ma(v30); ma(v31); ma(v32); ma(v33); ma(v4); ma(v5);
  se1 := root(d((1.0*n1*n1-n1-r)*(n1-1)*s*t,n1*r*r*r));
  se30 := (v30-n1)/root(tmm);
  se31 := (v31-n1)/root(tmm);
  se32 := (v32-n1)/root(tmm);
  se33 := (v33-n1)/root(tmm);
end; (* calcvhats1 *)

procedure calcvhats2(n1,m1: integer; s,r: real; v1,v3,v5: real;
          var vht,v6,v7,v8: real);
(* In the summation for the HT estimators,
   the global variables n and m are used. *)
var t,vvht,vv6,vv7,vv8,ca0,ca1,ca2,ca3: real;
    i,bb: integer;
begin
  t := r+s;
  ca1 := 1.0 - n1/v1;
  ca2 := 1.0 - n1/v3;
  ca3 := 1.0 - n1/v5;
  vvht := 0; vv6 := 0; vv7 := 0; vv8 := 0;
  for i:=n downto 1 do begin
    bb := inb[s0[i]];
    vvht := vvht + d(1.0, 1 - power(ca0,bb));
    vv6 := vv6 + d(1.0, 1 - power(ca1,bb));
    vv7 := vv7 + d(1.0, 1 - power(ca2,bb));
    vv8 := vv8 + d(1.0, 1 - power(ca3,bb));
   end;
  for i:=m downto 1 do begin
    bb := inb[s1[i]];
    vvht := vvht + d(1.0, 1 - power(ca0,bb));
    vv6 := vv6 + d(1.0, 1 - power(ca1,bb));
    vv7 := vv7 + d(1.0, 1 - power(ca2,bb));
    vv8 := vv8 + d(1.0, 1 - power(ca3,bb));
   end;
  if (vvht > n1+m1) then vht := vvht else vht := n1+m1;
  if (vv6 > n1+m1) then v6 := vv6 else v6 := n1+m1;
  if (vv7 > n1+m1) then v7 := vv7 else v7 := n1+m1;
  if (vv8 > n1+m1) then v8 := vv8 else v8 := n1+m1;
end; (* calcvhats *)


procedure snowstat;
(* calculates estimates and sampling variances ;
   vertices are indexed by i and j in population,
     and by k and kk in sample;
   uses global variables:
   n,m: sizes of initial sample and first wave sample;
   s0[k], s1[k]: vertex numbers for initial and first wave snowball sample;
   correspondence is i = s0[k], j = s1[kk];
   sample0: indicators in population of initial sample;
   sample1: indicators in population of first wave snowball sample;
   a[k] = out[i]: out-degrees;
   outball[i,kk]: j-numbers of vertices with an arc from vertex i,
      for kk = 1 to out[i]; *)

var i,j,k,kk,a0k,d0k,d1k,bb,nj,mj: integer;
    stat_s,stat_t,stat_r,stat_k,ca1,ca2,ca3: real;
    a0, b0, d0, d1, stat_kk: veci;
    vhat1,vhat3,vhat30,vhat31,vhat32,vhat33,vhat4,vhat5,
     vhatht,vhat6,vhat7,vhat8,ster1,ster3,ster30,ster31,ster32,ster33: real;
    ster1j,ster3j,ster30j,ster31j,ster32j,ster33j,ster4j,ster5j,
          ster6j,ster7j,ster8j: real;
    kj,rj,sj,v1j,v3j,v30j,v31j,v32j,v33j,v4j,v5j,
          vhtj,v6j,v7j,v8j,s1j,s30j,s31j,s32j,s33j: real;
    v1, v3, v30, v31, v32, v33, v4, v5, v6, v7, v8: vecrn;

procedure report(var ff: text);
begin
  writeln(ff);
  writeln(ff,'Statistics              ');
  writeln(ff,'n                 ', n:8);
  writeln(ff,'m                 ', m:8);
  writeln(ff,'k                 ', stat_k:8:0);
  writeln(ff,'s                 ', stat_s:8:0);
  writeln(ff,'r                 ', stat_r:8:0);
  writeln(ff);
  writeln(ff,'Estimates               ');
  writeln(ff,'v_1                 ', vhat1:8:1);
  writeln(ff,'v_3                 ', vhat3:8:1);
(*  writeln(ff,'v_3(3)              ', vhat33:8:1);  *)
  writeln(ff,'v_4                 ', vhat4:8:1);
  writeln(ff,'v_5                 ', vhat5:8:1);
if indegrees then begin
    writeln(ff,'v_6                 ', vhat6:8:1);
    writeln(ff,'v_7                 ', vhat7:8:1);
    writeln(ff,'v_8                 ', vhat8:8:1);
   end;
  writeln(ff);
  writeln(ff,'Standard errors (j = jackknife, m = modelbased);');
  writeln(ff, '(modelbased standard errors usually are of doubtful practical value).');
  writeln(ff,'se_1    m           ', ster1:8:1);
  writeln(ff,'se_1    j           ', ster1j:8:1);
  writeln(ff,'se_3    m           ', ster3:8:1);
  writeln(ff,'se_3    j           ', ster3j:8:1);
(*  writeln(ff,'se_3(3) m           ', ster33:8:1);
  writeln(ff,'se_3(3) j           ', ster33j:8:1); *)
  writeln(ff,'se_4    j           ', ster4j:8:1);
  writeln(ff,'se_5    j           ', ster5j:8:1);
  if indegrees then begin
    writeln(ff,'se_6    j           ', ster6j:8:1);
    writeln(ff,'se_7    j           ', ster7j:8:1);
    writeln(ff,'se_8    j           ', ster8j:8:1);
   end;
  writeln(ff);
  writeln(ff,'The conclusion section of Frank & Snijders (1994) favors estimates');
  writeln(ff,'v_3 and v_5, but notes that v_3 may be biased.');
  writeln(ff,'If v_3 and v_5 yield very different outcomes, v_3 is not reliable.');
  writeln(ff);
 end; (* report *)

begin
(* computation of several statistics *)
(* a0[k] is out-degree of k within initial sample;
   b0[j] is number of points in initial sample
   pointing to j *)
b0 := zeros;
for k:= 1 to n do begin
  a0k := 0;
  i := s0[k];
  for kk := a[k] downto 1 do begin
    j := outball[i,kk];
    if sample0[j] then Inc(a0k);
    Inc(b0[j]);
   end; (* kk *)
  a0[k] := a0k;
 end; (* k *)
if indegrees then begin
  for i:=1 to n do
    if (inb[i] < b0[i]) then begin
      writeln('There point more vertices in the initial sample to vertex ',i);
      writeln('than the given in-degree !');
      writeln('This is an inconsistency.');
      Waitkey;
     end;
  for i:=1 to m do
    if (inb[s1[i]] < b0[s1[i]]) then begin
   writeln('There point more vertices in the initial sample to vertex ',s1[i]);
   writeln('than the given in-degree !');
   Waitkey;
     end;
 end;
for k:=1 to n do Inc(inb[s0[k]]);
for k:=1 to m do Inc(inb[s1[k]]);

stat_r := 0; stat_t := 0; stat_k := 0;
stat_kk := zeros;
for k := 1 to n do begin
  i := s0[k];
  if (b0[i] >= 1) then begin stat_kk[k] := 1; stat_k := stat_k+1; end;
  stat_r := stat_r + a0[k];
  stat_t := stat_t + out[i];
  d0k := 0; d1k := 0;
  for kk := a[k] downto 1 do begin
    j := outball[i,kk];
    if (b0[j] <=1) then begin
      if (sample0[j]) then Inc(d0k) else Inc(d1k); end;
   end; (* kk *)
  d0[k] := d0k;
  d1[k] := d1k;
 end; (* k *)
stat_s := stat_t - stat_r;

(* computation of estimates *)
calcvhats1(n,m,stat_k,stat_s,stat_r,
    vhat1,vhat30,vhat31,vhat32,vhat33,vhat4,vhat5,
                 ster1,ster30,ster31,ster32,ster33);
vhat3 := v3_ex(n,m,stat_t,0,true);
if (stat_t - 0.1 > m) then ster3 := (vhat3-n)/root(stat_t - m)
   else ster3 := vhat3 - n;
if indegrees then
  calcvhats2(n,m,stat_s,stat_r,vhat1,vhat3,vhat5,vhatht,vhat6,vhat7,vhat8);

(* jackknifes *)
for k := 1 to n do begin
  i := s0[k];
  nj := n-1;
  mj := m + stat_kk[k] - d1[k];
  kj := stat_k - stat_kk[k] - d0[k];
  rj := stat_r - a0[k] - b0[i];
  sj := stat_s - a[k] + a0[k] + b0[i];
  calcvhats1(nj,mj,kj,sj,rj,
          v1j,v30j,v31j,v32j,v33j,v4j,v5j,s1j,s30j,s31j,s32j,s33j);
  v3j := v3_ex(nj,mj,sj+rj,vhat3,false);
  if indegrees then
    calcvhats2(nj,mj,sj,rj,v1j,v3j,v5j,vhtj,v6j,v7j,v8j);
  v1[k] := v1j;
  v3[k] := v3j;
(*  v30[k] := v30j;
  v31[k] := v31j;
  v32[k] := v32j; *)
  v33[k] := v33j;
  v4[k] := v4j;
  v5[k] := v5j;
  if indegrees then begin
(* Now subtract for the HT estimators
   the terms that should not have been added *)
    ca1 := 1.0 - nj/v1j;
    ca2 := 1.0 - nj/v3j;
    ca3 := 1.0 - nj/v5j;
    if (stat_kk[k] <= 0) then begin
      v6j := v6j - d(1.0, 1-power(ca1,inb[i]));
      v7j := v7j - d(1.0, 1-power(ca2,inb[i]));
      v8j := v8j - d(1.0, 1-power(ca3,inb[i]));
     end;
    for kk := out[i] downto 1 do
      if (b0[outball[i,kk]] <= 1) then begin
         bb := inb[outball[i,kk]];
         v6j := v6j - d(1.0, 1-power(ca1,bb));
         v7j := v7j - d(1.0, 1-power(ca2,bb));
         v8j := v8j - d(1.0, 1-power(ca3,bb));
        end;
    v6[k] := v6j;
    v7[k] := v7j;
    v8[k] := v8j;
   end; (* if indegrees *)
 end; (* k *)
jack(v1,ster1j);
jack(v3,ster3j);
(* jack(v30,ster30j);
jack(v31,ster31j);
jack(v32,ster32j); *)
jack(v33,ster33j);
jack(v4,ster4j);
jack(v5,ster5j);
if indegrees then begin
  jack(v6,ster6j);
  jack(v7,ster7j);
  jack(v8,ster8j);
 end;

report(output);
report(outf);

end; (* snowstat*)

begin (* main *)
initialise;
readfile;
snowstat;
close(outf);
waitkey;
end.
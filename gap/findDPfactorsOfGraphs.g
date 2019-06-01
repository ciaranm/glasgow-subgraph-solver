LoadPackage("digraph", false);

input:=InputTextUser();
firstLine:=ReadLine(input);
numOfNodes:=Int(SplitString(firstLine, "\n"," ")[1]);

digraphInput:=[];

for i in [1..numOfNodes] do
  line:=ReadLine(input);
  neighbourList:=List(SplitString(line, "", " \n"), x->Int(x));
  Remove(neighbourList,1);
  Add(digraphInput,neighbourList+1 );
od;

for x in digraphInput do
    IsHomogeneousList(x);
od;

gamma:=Digraph(digraphInput);

if not IsSymmetricDigraph(gamma) then
  gamma:= DigraphSymmetricClosure(gamma);
fi;

G:=AutomorphismGroup(gamma);
WriteAll(OutputTextUser(), Concatenation("--pattern-automorphism-group-size ", String(Size(G))));

for i in [1..numOfNodes] do
  orb := Orbit(Stabilizer(G,[1..i-1], OnTuples), i);
  for j in orb do
    if i <> j then
      WriteAll(OutputTextUser(), Concatenation(" --pattern-less-than '", String(i-1), "<", String(j-1), "'"));
    fi;
  od;
od;
WriteLine(OutputTextUser(), "");
QUIT_GAP();

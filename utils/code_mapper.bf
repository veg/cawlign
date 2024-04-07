LoadFunctionLibrary ("TemplateModels/chooseGeneticCode.def");

map = {64,1};

j = 0;
for (i,k; in; genetic_code.DefineCodonToAAMapping (_Genetic_Code)) {
   map[j] = k;
   j += 1;
}

fprintf (stdout, Join (",", map));
import json
from operator import itemgetter

def camel_case_split(str):
    words = [[str[0]]]
 
    for c in str[1:]:
        if c.isupper():
            words.append(list(c))
        else:
            words[-1].append(c)
 
    return [''.join(word) for word in words]

with open("benchmarking_query.json") as f:
    benchmarking_query=json.load(f)
benchmarking_query_by_sg = [[] for _ in range(230)]
for entry in benchmarking_query:
    benchmarking_query_by_sg[int(entry[0].split('_')[3])-1].append(entry)

models_required = set()
with open("test_generator.json","w") as f_tg, open("test_commands.sh","w") as f_tc:
    for i,sg_category in enumerate(benchmarking_query_by_sg):
        sg_category.sort(key=itemgetter(2))
        if len(sg_category)>0:
            model_name = sg_category[0][1]
            models_required.add(model_name)
            test_name = sg_category[0][0]
            test_name_split = test_name.split("_")
            prototype_label = "_".join(test_name_split[1:-5])
            species_string = test_name_split[-5]
            stoichiometric_species = camel_case_split(species_string)
            elastic_test_name_pattern = "ElasticConstants_"+prototype_label+"_"+species_string+"*"
            f_tc.write("pipeline-run-pair %s %s -v\n"%(elastic_test_name_pattern,model_name))
            f_tg.write(json.dumps({
                "stoichiometric_species": stoichiometric_species, 
                "prototype_label": prototype_label, 
                "crystal_genome_test_args": {},
            })+"\n")

with open("install_commands.sh","w") as f:
    for model in models_required:
        f.write("kimitems install -D %s\n"%model)

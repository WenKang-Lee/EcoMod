import random, math

#STEP 1: External & Internal Split; Create Presence/Absence Dict
## filter presence to retain uniques
presence_coord = [(35.2, 11.7), (-42.6, -32.1), (18.4, -6.3), (3.11, 8.12), (18.4, -6.3)]
''' Generate Pseudo absence, probably through incorporation of GIS. '''
absent_coord = [None]
uiq_coord = list(set(presence_coord))
''' here we perform conversion of coordinations into pixel index on the map for model assessment.'''
presence_idcs = [(4,7), (1,22), (9,6), (7,7)]
absent_idcs = [(12,6), (0,37), (9,9), (11,3)]
''' And we use these indices to search through and obtain the explanatory value.
Then we compund a tuple for each presence/ absence with the first element being the P/A bit.
Format example: (P/A_bit, env_val1, env_val2, env_val3)
! Maybe we should sort the env vars so that the first env element is the one starts with closest to first alphabet.
And we could record this sequence in a tuple. e.g.'''
env_idices = ("elev", "humi", "temp")
'''And while we're loop through the env raster files, we record the min n max of each env var for later use, as in:'''
minmaxes = ((3, 17), (11.1, 19.7), (6.43, 22.41))
'''Then we compound them in a list and shuffle the deck. Let's call it data_tuples. '''
data_tuples = [(1, 15.2, 22.1), (0, 31.6, 7.5), (1, 20.2, 11.5)]

## let's get some geostatistics from the env raster files beforehand such as:
env_1 = 17,1, 25.6
env_2 = 33.96, 26.12

## if we're doing 50/50 split
train_ds = random.sample(data_tuples, int(len(data_tuples)/2))
test_ds = [_ for _ in data_tuples if _ not in train_ds]

### sacred place
unq_i = 1; model_syn = {}; model_score = {}; tent_dict = {}

'''sub-functions repository'''
def var_roller(syn_holder, mm_i=None, mod_type=None):
    if mod_type == "rule":
        min_x, max_y = min(minmaxes[mm_i]), max(minmaxes[mm_i])
        range1 = random.uniform(min_x, max_y)
        roller1 = random.random()
        if roller1 < 0.33: # smaller <
            syn_holder.append(f"and psc_abs[{mm_i+1}] < {range1}") #! psc_abs
        elif roller1 < 0.66: # larger >
            syn_holder.append(f"and psc_abs[{mm_i+1}] > {range1}")
        else: # between a and b
            range2 = random.uniform(min_x, max_y)
            syn_holder.append(f"and {min(range1,range2)} < psc_abs[{mm_i+1}] < {max(range1, range2)}")
    
    elif mod_type == "log":
        BN = random.uniform(-3, 3)
        syn_holder.append(f"+ {BN}*psc_abs[{mm_i+1}]") #subsequent coefficients, BN
        
    else:
        raise KeyError("'rule' or 'log' only")
    
def substitution(mod, syn, tent_dict=None, insert=False, delete=False):
    global unq_i
    copy_syn = syn
    null_i = [i for i, x in enumerate(syn) if x == "and 1" or "+ 0"] #**
    var_i = [i for i in range(1,len(minmaxes)) if i not in null_i] # skip index 0 constant; mutual independent**
    
    if insert and delete:       #substitute checker
        if null_i and var_i:
            pass
        else:
            null_i = []; var_i = []
    
    if insert:
        var_clip = []
        if null_i:
            new_vi = random.choice(null_i)-1
            subt_type = "rule" if mod[0] == "R" else "log"
            var_roller(syn_holder=var_clip, mm_i=new_vi, mod_type=subt_type)
            copy_syn[new_vi] = var_clip[0]
        
    if delete:
        if var_i:
            getrid_vi = random.choice(var_i)
            copy_syn[getrid_vi] = "and 1" if mod[0] == "R" else "+ 0"
        
    tent_dict[f"{mod[0]}_{unq_i}"] = copy_syn
    unq_i += 1

def dict_trim():
    for selected_key in set(model_syn).intersection(model_score):
            if selected_key not in model_syn:
                del model_syn[f"{selected_key}"]

#STEP 2: Initiates pool w/ 4 types of rules
## if we're doing equal size types)
def model_factory(evo=False, n=50):
    global model_syn, model_score, unq_i, tent_dict
    for _ in range(n):
        mod_type = "rule" if random.random() <= 0.5 else "log"
        words = []; words.append("1") if mod_type == "rule" else words.append(str(random.randint(-20, 20)))
        for i in range(len(minmaxes)):
            # include corresponding env var
            var_roller(words, mm_i=i) if random.random() > 0.5 else words.append("and 1") if mod_type == "rule" else words.append("+ 0")
        prefix = "R" if mod_type == "rule" else "L"
        (model_syn if not evo else tent_dict)[f"{prefix}_{unq_i}"] = words

def model_evaluate(iter_mod=None):
    for mod, syn in iter_mod.items():
        model_score[mod] = (0,0,0,0,None,None,None) #TP, FP, FN, TN, Accuracy, CE, OE
        model_x = " ".join(syn)
        for psc_abs in train_ds:    # loop through every single grid
            if mod[0] == "L": #if it's Log. Regression
                logit = eval(model_x)
                prob = 1/(1+math.exp(-logit))
                outcome = 1 if prob > 0.5 else 0 #let's say we use 50% as minimum presence mark
            else:
                outcome = 1 if eval(model_x) else 0
            # Evaluation
            if outcome == 1 and psc_abs[0] == 1: #True Positive
                model_score[mod][0] += 1
            elif outcome > psc_abs[0]: #False Positive
                model_score[mod][1] += 1
            elif outcome < psc_abs[0]: #False Negative
                model_score[mod][2] += 1
            elif outcome == 0 and psc_abs[0] == 0: #True Negative
                model_score[mod][3] += 1

        #fitness function
            total_sum = model_score[mod][0] + model_score[mod][1] + model_score[mod][2] + model_score[mod][3]
            accuracy = (model_score[mod][0] + model_score[mod][3])/total_sum
            model_score[mod][-3] = accuracy
            comm_err = model_score[mod][0]/(model_score[mod][0] + model_score[mod][1])
            model_score[mod][-2] = comm_err
            om_err = model_score[mod][0]/(model_score[mod][0] + model_score[mod][2])
            model_score[mod][-1] = om_err

def nature_select(CET=0.3, OET=0.1):  #commission and omission errors limit
    global unq_i
    generation = 1
    increment = 1                   #to initiate the loop
    while increment >= 0.1 and generation < 101:
        model_score = {mod: score for index, (mod, score) in enumerate(sorted(model_score.items(), key=lambda x: x[1][-3])) if index <= round(len(model_score)/2)}
        ances_apex = next(iter(model_score.items()))[1][-3]
        #^^^ retain first-half in accuracy
        dict_trim()
        
        # mutation, crossover, and integration
        for _, (mod, syn) in enumerate(model_syn):
            dice = random.randint(1,6)
            if dice == 1:                           # crossover
                parent_b = "to_initiate_the_crossover"
                parent_a = mod
                while parent_a == parent_b or parent_a[0] != parent_b[0]:
                    parent_b = random.choice(list(model_syn))
                crossing_point = random.randint(0, len(minmaxes))
                tent_dict[f"{parent_a[0]}_{unq_i}"] = model_syn[parent_a][:crossing_point] + model_syn[parent_b][crossing_point:]
                unq_i += 1
                tent_dict[f"{parent_a[0]}_{unq_i}"] = model_syn[parent_b][:crossing_point] + model_syn[parent_a][crossing_point:]
                unq_i += 1
                
            elif dice == 2:                         # insertion
                substitution(mod, syn, tent_dict=tent_dict, insert=True)
                
            elif dice == 3:                         # deletion
                substitution(mod, syn, tent_dict=tent_dict, delete=True)
                
            elif dice == 4:                         # substitution
                substitution(mod, syn, tent_dict=tent_dict, insert=True, delete=True)
                
            elif dice == 5:                         # integration (OR)
                parent_b = "kick start model merger"
                parent_a = mod
                while parent_a == parent_b:
                    parent_b = random.choice(list(model_syn))
                tent_dict[f"H_{unq_i}"] = model_syn[parent_a] + ["or"] + model_syn[parent_b] # H for hybrid. end species in evolutionary race.
                unq_i += 1
                
            else: #dice == 6                        # creation
                model_factory(evo=True, n=1)
            
        model_syn.update(tent_dict)
        model_evaluate(iter_mod=tent_dict)
        increment = next(iter(model_score.items()))[1][-3] - ances_apex
        tent_dict = {}
        generation += 1
    
    model_score = {mod: score for mod, score in sorted(model_score.items(), key=lambda x: x[1][-3])}
    checker_i = 0
    while checker_i < len(model_score) and checker_i < 10:
        del_key = list(model_score)[checker_i]
        if model_score[f"{del_key}"][-2] > CET or model_score[f"{del_key}"][-1] > OET:
            del model_score[f"{del_key}"]
        checker_i += 1
    model_score = {mod: score for i, (mod, score) in model_score.items() if i < 10}
    dict_trim()
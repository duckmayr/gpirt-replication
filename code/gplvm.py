import numpy as np
import GPy

from GPy.models.bayesian_gplvm_minibatch import BayesianGPLVMMiniBatch

full_responses = np.genfromtxt('data/npi-responses.csv', delimiter=',')

seeds = [969167, 188942, 134058, 124022, 685285, 226318, 365209, 648795,
         985797, 193627, 569692, 589449, 832867, 497690, 402858, 583422,
         183204, 883281, 669543, 277324]

print("Running GPLVM replicates")
for iteration in range(20):
    print("\tRunning replicate " + str(iteration+1) + "\n")
    np.random.seed(seeds[iteration])
    altered_responses = np.copy(full_responses)
    remove_idx = np.random.choice(np.arange(altered_responses.size),
                                  replace=False,
                                  size=int(altered_responses.size * 0.2))
    outname = 'model-output/GPLVM-iter-' + str(iteration) + '-idx.csv'
    np.savetxt(outname, remove_idx, delimiter = ',', newline = '\n')
    altered_responses[np.unravel_index(remove_idx, altered_responses.shape)] = np.nan
    m = BayesianGPLVMMiniBatch(altered_responses, 1, missing_data=True)
    m.optimize(messages=0, max_iters=5e3)
    pred = m.predict(full_responses)
    outname = 'model-output/GPLVM-iter-' + str(iteration) + '-mean.csv'
    np.savetxt(outname, pred[0], delimiter = ',', newline = '\n')
    outname = 'model-output/GPLVM-iter-' + str(iteration) + '-var.csv'
    np.savetxt(outname, pred[1], delimiter = ',', newline = '\n')
print("")


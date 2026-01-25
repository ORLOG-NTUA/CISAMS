import numpy as np
import pyjags
import xarray as xr


##################### Bayesian BWM main function
def BayesianBWM(ab,aw):

  #kostas edit: 8.8.2025
  #sett random seed for replicability
  ###import numpy as np
  ###np.random.seed(1977)
  #end kostas edit: 8.8.2025

  dm_no, c_no = aw.shape

  wStar_prior = (1/10000000) * np.ones((c_no,));

  data_list = {
      "AW": aw,
      "AB": ab,
      "wStarPrior":wStar_prior,
      'n':np.sum(ab,axis=1),
      'm':np.sum(aw, axis=1),
      'k':dm_no,
      'c':c_no,
  }

  BayesianBWM_model = '''
    model {
      for(i in 1:k){
        AB[i,1:c] ~ dmulti(1/w[i,1:c],n[i])
        AW[i,1:c] ~ dmulti(w[i,1:c],m[i])
        w[i,1:c] ~ ddirch(beta*wStar)
      }

        beta ~ dgamma(0.01 , 0.01)
        wStar ~ ddirch(beta2*theta)
        theta ~ ddirch(wStarPrior)
        beta2 ~ dgamma(0.01,0.01)

    }
    '''
  parameters = ['w', 'wStar']
  variables = parameters + ['log_like']


  chain_no = 4
  sample_per_chain = 5000
  sample_burnin = 500

  jags_BayesianBWM = pyjags.Model(
      code=BayesianBWM_model,
      data=data_list,
      chains=chain_no,
      threads=4,
      chains_per_thread=1
  )

  samples = jags_BayesianBWM.sample(sample_per_chain + sample_burnin, vars=parameters)
  samples = pyjags.discard_burn_in_samples(samples, burn_in=sample_burnin)

  wStar_samples = np.reshape(samples['wStar'], [c_no, chain_no*sample_per_chain])
  wStar_mean = np.mean(wStar_samples, axis=1)
  print("The average of the aggregated weight distribution is:", wStar_mean)

  return wStar_mean, wStar_samples.T


##################### Visualization of credal ranking
def CredalRanking(wStar, criteria_name):
  import math
  import networkx as nx
  import matplotlib.pyplot as plt
  from matplotlib.colors import Normalize
  import seaborn as sns
  import matplotlib.cm as cmx


  def roundUp(n, d=2):
      d = int('1' + ('0' * d))
      return math.ceil(n * d) / d

  # np.mean(w_star,axis=0)
  n_c = wStar.shape[1]
  w = np.zeros((n_c, n_c))

  for i in range(n_c):
    for j in range(i+1, n_c):
      w_ij = np.sum(wStar[:,i] > wStar[:,j])
      w_ji = np.sum(wStar[:,i] < wStar[:,j])
      if w_ij > w_ji:
        w[i,j] = roundUp((w_ij / wStar.shape[0]), 2)
      else:
        w[j,i] = roundUp((w_ji / wStar.shape[0]), 2)


  colors = sns.color_palette(None, n_c)

  plt.figure(figsize=(12,12))

  # Create DiGraph from w
  G = nx.convert_matrix.from_numpy_array(w, create_using=nx.DiGraph)

  # Use spring_layout to handle positioning of graph
  layout = nx.kamada_kawai_layout(G)


  # Color mapping
  colors = plt.cm.jet #sns.color_palette(None, n_c)
  cNorm  = Normalize(vmin=0, vmax=n_c-1)
  scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=colors)
  values = [node for node in G.nodes()]

  # Using a figure to use it as a parameter when calling nx.draw_networkx
  f = plt.figure(1)
  ax = f.add_subplot(1,1,1)
  for label in criteria_name:
      ax.plot([0],[0],color=scalarMap.to_rgba(criteria_name[label]),label=label)

  # Draw the graph using the layout
  nx.draw(G, layout, with_labels=True, node_size=2000, edgecolors='black', node_color=values,ax=ax)

  # Get weights of each edge and assign to labels
  labels = nx.get_edge_attributes(G, "weight")

  # Draw edge labels using layout and list of labels
  nx.draw_networkx_edge_labels(G, pos=layout, edge_labels=labels)

  # Setting it to how it was looking before.
  plt.axis('off')
  f.set_facecolor('w')

  ###plt.legend()

  f.tight_layout()
  ###plt.show()
  plt.savefig('Circularity_1f.jpeg')
  return w


#A_B =  np.array([
#       [1, 7, 3, 4, 2, 6, 9, 5, 3, 4, 5, 4],
#       [3, 2, 1, 3, 4, 5, 6, 6, 7, 7, 8, 9],
#       [6, 5, 6, 5, 6, 6, 7, 7, 7, 8, 7, 6],
#       [3, 3, 5, 1, 2, 2, 4, 2, 3, 4, 2, 2],
#       [2, 1, 3, 4, 2, 4, 2, 6, 2, 6, 4, 4],
#       [4, 1, 3, 5, 6, 7, 3, 4, 8, 4, 2, 3],
#       [3, 3, 2, 3, 1, 3, 6, 4, 2, 3, 3, 3],
#       [1, 4, 5, 2, 7, 8, 6, 3, 9, 8, 2, 4]])
A_B = np.loadtxt('data_CIRCULARITY_Market_A_B.txt')
check1 = np.loadtxt('check_CIRCULARITY_Market_rev.txt')
A_B = A_B[check1==1,:]
print(A_B)

#A_W =  np.array([
#      [ 9, 2, 7, 6, 7, 5, 1, 6, 5, 5, 6, 4],
#      [ 9, 8, 7, 6, 5, 5, 4, 4, 3, 3, 1, 2],
#      [ 8, 8, 7, 8, 7, 7, 8, 8, 8, 7, 7, 8],
#      [ 2, 3, 1, 5, 4, 5, 3, 4, 2, 2, 3, 5],
#      [ 5, 6, 3, 2, 5, 2, 5, 1, 5, 1, 4, 4],
#      [ 5, 8, 6, 4, 3, 2, 6, 5, 1, 5, 7, 6],
#      [ 4, 4, 5, 4, 6, 4, 1, 3, 5, 4, 4, 4],
#      [ 9, 5, 4, 7, 3, 2, 6, 8, 1, 3, 7, 6]])
A_W = np.loadtxt('data_CIRCULARITY_Market_A_W.txt')
check1 = np.loadtxt('check_CIRCULARITY_Market_rev.txt')
A_W = A_W[check1==1,:]
print(A_W)

criteria_name = {'C1': 0, 'C2':1, 'C3':2, 'C4':3,'C5':4, 'C6':5, 'C7':6, 'C8':7, 'C9':8}
#criteria_name = {'C1': 1, 'C2':2, 'C3':3, 'C4':4,'C5':5, 'C6':6, 'C7':7, 'C8':8, 'C9':9}

# kostas edit: 8.8.2025
# sett random seed for replicability
###import numpy as np

###np.random.seed(1977)
# end kostas edit: 8.8.2025

wStar_mean, wStar_samples = BayesianBWM(A_B, A_W)

probabilities = CredalRanking(wStar_samples, criteria_name)

print(probabilities)

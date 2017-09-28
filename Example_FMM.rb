# Example RevBayes script to set up a finite mixture model (FMM) for site-specific rates of evolution

# --> Setup <--

# Read in data (amino acid, in this case)
data <- readDiscreteCharacterData("myData.nex")

# Read in tree
# Note: This version of the analysis uses a fixed topology to improve mixing for site rates.
myTree <- readBranchLengthTrees("myTree.tre")

# --> Get some useful variables from the data. We need these later on. <--

n_species <- data.ntaxa()		# Number of species
n_sites <- data.nchar()			# Number of sites
taxa <- data.taxa()				# Taxon Names
n_branches <- 2*n_species-3		# Number of branches in (unrooted) tree


# --> Set Move and Monitor Indices <--

mvi = 0
mni = 0

# --> Model <--

# All rate categories evolve according to the same transition matrix
Q := fnJones()		# JTT matrix for protein evolution

NUM_MIX_CATS <- 6	# define the number of mixture elements

# For each mixture element, we need to create a stochastic variable corresponding
# to the tree length for all sites associated with that element. 

for (i in 1:NUM_MIX_CATS) {

	#Gamma prior on TL, in line with Rannala et al (2011)
    site_rate_value[i] ~ dnGamma(1,0.1)	
    
    # Add a scaling move for each tree length
    moves[++mvi] = mvScale(site_rate_value[i],lambda=1.0,tune=true,weight=2)
}

# Define a Dirichlet hyperprior on the probabilities of assigning sites to mixture components
mixture_probs ~ dnDirichlet(rep(1,NUM_MIX_CATS))

# Adding moves for assignment probabilities
moves[++mvi] = mvSimplexElementScale(mixture_probs)
moves[++mvi] = mvSimplex(mixture_probs,alpha=0.1,weight=1.0)

# For each site, create a stochastic variable that draws from the mixture.

for (i in 1:n_sites) {

	# Create and initialize stochastic draw from mixture for site i.
    siteRates[i] ~ dnMixture(site_rate_value,mixture_probs)

	# Adding moves for site i's mapping variable
    moves[++mvi] = mvMixtureAllocation(siteRates[i], weight=1.0)
    moves[++mvi] = mvGibbsMixtureAllocation(siteRates[i], weight=1.0)
	# Note: Because we are setting up moves for each site independently, the vast majority
	#       of MCMC proposals will correspond to these.
}

# Create a fixed variable for the unrooted tree topology
psi <- myTree[1]

# Create stochastic variables for the branch lengths on the tree
br_lens ~ dnDirichlet(rep(1,n_branches))

# Adding moves for the branch lengths
moves[++mvi] = mvSimplexElementScale(br_lens,alpha=0.1,weight=1.0)
moves[++mvi] = mvSimplex(br_lens,alpha=0.1,weight=1.0)

# Define a deterministic variable for the tree length (sum of branch lengths)
TL := sum(br_lens)

# Define a deterministic variable for the overall phylogeny (topology + branch lengths)
phylogeny := treeAssembly(psi, br_lens)

for (i in 1:n_sites) {

	# This loop is slow, so this is a progress indicator.
	if (i % 10 == 0) { i }

	# Defining a continuous-time Markov chain for each site, based on the tree length
	# for that site drawn from the mixture.
    phyloSeq[i] ~ dnPhyloCTMC(tree=phylogeny, Q=Q, branchRates=siteRates[i], type="AA")

	# To clamp the site to the appropriate CTMC, we first exclude all data...
    data.excludeAll()

	# ...then re-include only that site.
    data.includeCharacter(i)

	# Now we clamp the CTMC
    phyloSeq[i].clamp(data)
}

# We create a model object and assign it a stochastic node from our graphical model
mymodel = model(phylogeny)

# We can create various types of monitors to record analysis progress
monitors[++mni] = mnModel(filename="myAnalysis.log",printgen=1, separator = TAB)
monitors[++mni] = mnFile(filename="myAnalysis.trees",printgen=1, separator = TAB, psi)
monitors[++mni] = mnScreen(printgen=1, TL)

# We create our MCMC object and then start our analysis
mymcmc = mcmc(mymodel, monitors, moves)
mymcmc.run(generations=100000)


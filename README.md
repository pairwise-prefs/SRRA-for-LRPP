#The what?
Implementing the SRRA method for LRPP.
Including a benchmark data.

- - -
#code details
Implementated in Ruby (http://www.ruby-lang.org/en/) version 1.9.2
Dependecies: 

##Preliminaries
 Gem dependencies: rmagick, tsort  -- note these dependenceis are used only for pretty visualize/analyse data; feel free to remove them
 Native c implementations: under the ext/ folder; for each of its subfolders perform: ruby extconf.rb; then "make". copy the libsvm result to lib/  

#Usage example
 benchmark = JSON.parse File.new( DATA_XXX_JSON,"r").readlines.join( " " )
 
 benchmark_tag = benchmark["tag"]
 
 problem = benchmark["datas"][ PROBLEM_INDEX ]
 
 pref = EpsSmooth::PreferenceOracle.new problem["W"]
 
 learner = EpsSmooth::SVMLearner.new problem["X"]
 
 m = 100 

 agg_sampled = [ ]

 pivot= learner.rANDOM( (0..xs.size-1).to_a, EpsSmooth::Perm.new((0..xs.size-1).to_a.shuffle), m )

 5.times do |i|
   
   pivot_cost = pref.kendalTauish_cost( EpsSmooth::Perm.new( svmranker.ranking_of( xs ) ) )
   
   puts "iter #{i} cost #{pivot_cost}"

   pivot,agg_sampled = learner.improve_all_data( EpsSmooth::Perm.new( pivot.ranking_of( xs ) ), pref, m, agg_sampled, beta=1.0 ); 
 end



require './libsvm.bundle'


def format( u )
  u.inject([ [ ], 0]){ |aux, v| node = Libsvm::Node.new; node.index = aux.last; node.value = v; aux[0] << node; aux[1] += 1; aux }.first
end

# This library is namespaced.
problem = Libsvm::Problem.new
parameter = Libsvm::SvmParameter.new

parameter.cache_size = 1 # in megabytes

parameter.eps = 0.001
parameter.c = 10

examples = [ [1,0,1], [-1,0,-1] ].map {|ary| format(ary) }
labels   = [1, -1]
weights  = [1, 1]
problem.set_examples(labels, examples, weights)

model = Libsvm::Model.train(problem, parameter)



# tests #
[ [0.5,0,0.5], [-0.5,0,-0.5], [0,0,0], [100,1,10] ].each do |x|
  test_x = format x

  pred = model.predict( test_x )

  soft_pred = model.soft_predict( test_x )

  svs_info = model.get_svs_info
  
  puts "svs: #{svs_info}"
  
  puts "Example #{x} - Predicted #{pred} , soft #{soft_pred}"
end



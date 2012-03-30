require File.dirname(__FILE__)+'/libsvm.bundle'

module SVMRanker

class Svm
  def initialize(data)
    data[:instances].map! { |u| format(u) }
    @model = svm_learn data
  end

  def pref(u, v)
    @model.predict format( vect_diff( u,v ) )
  end

  def soft_pref(u, v)
    ( @model.soft_predict format( vect_diff( u,v ) ) ).first
  end

  def ranking_of( alternatives )
   degs = (0..alternatives.size-1).inject( {} ) do |deg, u|
      (u+1..alternatives.size-1).inject( deg ) do |adeg, v|
        uv_pref = pref( alternatives[u], alternatives[v] )
        adeg[u] = (adeg.key?( u ) ? adeg[u] : 0) +  ( uv_pref < 1.0 ? 0 : 1 )
        adeg[v] = (adeg.key?( v ) ? adeg[v] : 0) +  ( uv_pref > 0.0 ? 0 : 1 )
        adeg
      end
    end

    degs.sort { |u,v| v.last <=> u.last }.map { |kv| kv.first }
  end

  def hinge_loss( u, v, pref_u_v )
    [0, 1 - (2*pref_u_v - 1) * soft_pref( u,v )].max
  end

#  :private
#  was:  (u.keys - v.keys + v.keys).sort.inject({}) { |uv, fi| uv[fi.to_i] = (u[fi].nil? ? 0.0 : u[fi]) - (v[fi].nil? ? 0.0 : v[fi]); uv }
#  current implementation is about 2 times faster, and much uglier -- but what the heck..
  def vect_diff2(u,v)
    ukeys = u.keys
    vkeys = v.keys
    inx_u = 0
    inx_v = 0
    diff = {}
    while inx_u<ukeys.size || inx_v<vkeys.size
      ufi = inx_u<ukeys.size ? ukeys[ inx_u ] : vkeys[ inx_v ] + 10
      vfi = inx_v<ukeys.size ? vkeys[ inx_v ] : ukeys[ inx_u ] + 10
      if ufi < vfi
        diff[ ufi.to_i ] = u[ ufi ]
        inx_u += 1
      else if ufi > vfi
             diff[ vfi.to_i ] = -v[ vfi]
             inx_v += 1
           else # ufi == vfi
             diff[ ufi.to_i ] = u[ufi] - v[ufi] 
             inx_v += 1
             inx_u += 1
           end
      end
    end
    diff
  end

  def vect_diff(u,v)
    ukeys = u.keys
    vkeys = v.keys
    inx_u = 0
    inx_v = 0
    diff = {}
    while inx_u<ukeys.size || inx_v<vkeys.size      
      ufi = inx_u<ukeys.size ? ukeys[ inx_u ] : nil
      vfi = inx_v<vkeys.size ? vkeys[ inx_v ] : nil

      fi = [ufi,vfi].compact.min{|a,b| a.to_i <=> b.to_i}

      diff[ fi.to_i ] = (u[ fi ]||0.0) - (v[ fi ]||0.0)
      inx_u += 1 if ufi == fi
      inx_v += 1 if vfi == fi
    end
    diff
  end




  def format( u )
    u.map {|f| node = Libsvm::Node.new; node.index = f.first.to_i; node.value = f.last; node}
  end

  def svm_conf
    parameter = Libsvm::SvmParameter.new

    parameter.cache_size = 256 # in megabytes
    parameter.kernel_type = 0 # linear kernel
    parameter.gamma = 0.001
    parameter.eps = 0.01
    parameter.c = 1
#    parameter.h = 1
    parameter
  end

  def svm_learn( svm_data )
    problem = Libsvm::Problem.new
    problem.set_examples svm_data[:labels], svm_data[:instances], svm_data[:weights]
    Libsvm::Model.train problem, svm_conf
  end
end

end # module

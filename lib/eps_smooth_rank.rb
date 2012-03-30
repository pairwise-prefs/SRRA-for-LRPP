require File.dirname(__FILE__)+'/svm_ranker.rb'
require 'set'
require 'rmagick'
require 'tsort'

require File.dirname(__FILE__)+'/../ext/eps_srra_native.bundle'
require File.dirname(__FILE__)+'/../ext/dist_aux.bundle'

class Hash
  include TSort
  alias tsort_each_node each_key
  def tsort_each_child(node, &block)
    fetch(node).each(&block)
  end
end

module EpsSmooth
# --
class Perm
 
  def initialize( order )  
    raise "Not a valid order (repitition): #{order}" unless valid? order
    @pi = order    
    @inv_pi = ( @pi.inject( {:aux =>(0..@pi.size-1).to_a, :res => Array.new(@pi.size) }) { |zz, r| zz[:res][r] = @pi.index(r); zz } ) [:res]    
  end

  def rank(alternative_inx)
    @pi[ alternative_inx ]
  end

  def inv_rank(rank)
    @inv_pi[ rank ]
  end

  def pref(u,v)
    rank(u) < rank(v) ? 1 : 0 
  end
  
  def size
    @pi.size
  end

  def move(u, new_inx)
    inx = rank u
    aux = @pi.map do |vi|  
                val = vi
                val = vi - 1 if inx < new_inx && vi > inx && vi <= new_inx 
                val = vi + 1 if inx > new_inx && vi >= new_inx && vi < inx 
                val
    end 
    aux[u] = new_inx
    Perm.new aux
  end

  :private 
  def valid?( perm_arr )
    Set.new( perm_arr ).size == perm_arr.size 
  end

end

# --
class Sampler
  def self.sample( pi, quota )
    u_quota = quota / pi.size
    u_dists = EpsSrraNative::rank_pi_distributions( pi )   
    agg_samples = []
    u_dists.each_with_index do |d_u, u|
      agg_samples += EpsSrraNative::rank_sample_for_u( u, d_u, pi, u_quota )      
    end
    agg_samples
  end

  def self.sample2( pi, quota, beta=1.0 )
    u_quota = quota / pi.size
    u_dists = EpsSrraNative::rank_pi_beta_distributions( pi, beta )   
    
    agg_samples = []
    u_dists.each_with_index do |d_u, u|
      agg_samples += EpsSrraNative::rank_sample_for_u( u, d_u, pi, u_quota )      
    end
#    puts "sample2 -- #{agg_samples}"
#    puts ""
    agg_samples
  end

  def self.random_sample( vs, quota )
    u_quota = quota / vs.size
    z = vs.inject( [ ] ) do |agg_samples,u|       
      vs_no_u = vs - [u] 
      size_vs_no_u = vs_no_u.size
      u_quota.times { v = rand( size_vs_no_u ); v += 1 if v >= u; agg_samples << [ u, v, nil, 1.0/size_vs_no_u, size_vs_no_u.to_f/u_quota ] } #note we correct v's index 
      agg_samples
    end
    #debug "rand -- #{z}"
    #debug ""
    z
  end

  def self.rank_pi_beta_distributions( pi, beta )
    all_dists = []
    norm_u = []
    (0..pi.size-1).each do |u|
      norm_u[u] = 0
      (0..pi.size-1).each do |v|
        pr_uv = v==u ? 0.0 : 1.0 / ( (pi.inv_rank( 0 )-pi.inv_rank(22)).abs )**beta
      end
    end

  end


end
#--


# --
class Estimator 
  
  def initialize( pi, eps, hp_r=100 )
    @n = pi.size
    @coin_n = 1.0 #hp_c * Math.log2(@n) 
    @pi = pi
    @s = (0..hp_r-1).inject( [] ) do |agg, i|
      agg += do_sample
    end    
  end

  def the_sample
    @s
  end

  def sample_the_sample( m )
    @s.sample( m )
  end

  def estimation(sigma, h)
    ecost(sigma, h) - ecost(@pi, h)
  end

  def to_s
    "EpsSmooth::Estimator{ n=#{@n}, |sample|=#{@s.size}, coin-const=#{@coin_n}}"
  end

  :private
  def do_sample( )    
    (0..@n - 1).inject( [ ] ) do |aux, i|
      (i+1..@n - 1).inject(aux) do |agg_sample, j|
        agg_sample << { :pair => [ i, j ], :pref => (@pi.rank(i) < @pi.rank(j) ? 1 : 0), :coin => coin(i,j) }  if draw?(i, j) 
        agg_sample
      end
    end
  end  

  def ecost(sigma, h)
    @s.inject( 0.0 ) do |agg, x| 
      u=x[:pair].first
      v=x[:pair].last 
      agg + ( sigma.pref(u,v) * h[v][u] + sigma.pref(v,u) * h[u][v] ) * (1.0 / x[:coin])
    end
  end
  
  def draw?(i,j)
    rand < coin( i, j )
  end

  def coin(i, j)
    [1.0, @coin_n / ( @pi.inv_rank(i) - @pi.inv_rank(j) ).abs].min 
  end

  def debug
    "Estimator debug: \n#{to_s}\n\nsample=#{@s}\n"
  end
end

# --
class PreferenceOracle
  attr_reader :n, :query_cost

  def initialize( pref_W )
    @w = pref_W
    @n = @w.first.size
    @query_cost = 0
  end

  def query(u,v)
    @query_cost += 1
    @w[u][v]
  end  


  #.......
  def kendalTauish_cost( pi )
    (0..@n-1).inject( 0 ) { |agg, u| (u+1..@n-1).inject( agg ) { |agg, v| agg += 1 unless pi.pref(u,v) == @w[u][v]; agg } }      
  end

  def hinge_loss( svmranker, xs )
    (0..@n-1).inject( 0.0 ) { |agg, u| (u+1..@n-1).inject( agg ) { |aux, v| aux + svmranker.hinge_loss(xs[u],xs[v], @w[u][v]) + svmranker.hinge_loss(xs[v],xs[u],@w[v][u]) } }
  end

  def opt_kt_cost
    0.5 * (0..@w.size-1).to_a.inject( 0 ) { |agg, u|  (0..@w.size-1).to_a.inject( agg ) { |agg2,v| agg2+ (@w[u][v].nil? ? 0 : ( u < v ? (1-@w[u][v]).abs : (1-@w[v][u]).abs ))  }  }
  end  
  #........
  # Native versions:    -- benchmark shows it does not worth the effort..
  def native_kendalTauish_cost( pi )
    EpsSrraNative::rank_aggergate_upperW(@n){ |u,v| pi.pref(u,v) == @w[u][v] ? 0 : 1 }    
  end

  def native_opt_kt_cost
    0.5 * EpsSrraNative::rank_aggergate_upperW(@n){  |u,v| @w[u][v].nil? ? 0 : ( u < v ? (1-@w[u][v]).abs : (1-@w[v][u]).abs )  }  
  end    

  def native_hinge_loss( svmranker, xs )
    EpsSrraNative::rank_aggergate_upperW(@n){ |u,v| svmranker.hinge_loss(xs[u],xs[v], @w[u][v]) + svmranker.hinge_loss(xs[v],xs[u],@w[v][u]) }
  end
  #...........
  def strongly_connected_components
    aux = @w.inject( {:i => 0, :g => {} } ){ |a,r| a[:g][a[:i]] = []; r.each_with_index{|v,i| a[:g][a[:i]] << i if v == 1 }; a[:i] = a[:i]+1; a }
    aux[:g].strongly_connected_components
  end

  ONE_COLOR = "0,26,51"
  ZERO_COLOR = "122,149,255"
  def to_img( fname,out_dir='./',perm=Perm.new( (0..@n-1).to_a ) )
    img = Magick::Image.new(@w.first.size, @w.size)
    
    (0..perm.size-1).each do |u|
      row_index = perm.rank u
      (0..perm.size-1).each do |v|        
        column_index = perm.rank v
        img.pixel_color(column_index,row_index, "rgb(224,224,224)") if row_index == column_index
        img.pixel_color(column_index,row_index, "rgb(#{@w[u][v] < 1 ? ZERO_COLOR : ONE_COLOR})") unless row_index == column_index
      end
    end

    img.write("#{out_dir}#{fname}")
  end
end

# --
class SVMLearner
  def initialize(alternatives)
    @xs = alternatives
  end

  def improve_new(perm, oracle, m)
    build_svmrank( oracle, Sampler.sample(perm,m) )
  end

  def improve_all_data(perm, oracle, m, agg_sampled, beta=1.0)
    the_sample = Sampler.sample2(perm,m,beta) + agg_sampled
    ranker = build_svmrank( oracle, the_sample )
    [ ranker, the_sample ]
  end

  def combi_improve(pi, oracle, m, agg_sampled)
    the_sample = Sampler.sample2(pi,m) + agg_sampled
    w_edges = the_sample.inject({}){ |m,x| uv = x[0]<x[1] ? [ x[0],x[1] ] : [x[1],x[0]]; m[uv] = [ oracle.query(uv.first, uv.last), 0.0 ] unless m.has_key? uv; m[uv][1] = m[uv][1]+x.last; m } # hash: key=[u,v] (u<v); value=[ true-pref, agg-weight ]
    [sparse_local_improve(pi, w_edges), the_sample]
  end

  def combi_rANDOM_improve(pi, oracle, m)    
    the_sample = Sampler.random_sample( (0..pi.size-1).to_a, m )
    w_edges = the_sample.inject({}){ |m,x| uv = x[0]<x[1] ? [ x[0],x[1] ] : [x[1],x[0]]; m[uv] = [ oracle.query(uv.first, uv.last), 0.0 ] unless m.has_key? uv; m[uv][1] = m[uv][1]+x.last; m } # hash: key=[u,v] (u<v); value=[ true-pref, agg-weight ]
    sparse_local_improve(pi, w_edges)
  end

  def rANDOM( vs, oracle, m )
    build_svmrank( oracle, Sampler.random_sample( vs, m ) )
  end
  
  def build_svmrank(oracle, sample_with_dups)
    s_no_dups = sample_with_dups.inject({}){ |m,x| m[x] = 0 unless m.has_key? x; m[x] = m[x]+x.last; m } # hash: key=result-info; value=aggregated-weight       

    est = s_no_dups.keys.inject( { :labels =>[], :instances =>[], :weights =>[]} ) do |data,x|
      if rand > 0.5
        u = x[0]
        v =  x[1]
      else
        u = x[1]
        v =  x[0]
      end
      data[:instances] << vect_diff_fast_and_ugly( @xs[ u ], @xs[ v ] )
      data[:labels] << oracle.query( u,v )
      data[:weights] << s_no_dups[ x ]
      data
    end
    
    svmrank = SVMRanker::Svm.new est
  end

  
  def to_s
    "EpsSmooth::Learner"
  end

#  :private
  def sample_weight_kt(pi, edges)
    n = pi.size
    (0..n-1).inject( 0 ) { |agg, u| (u+1..n-1).inject( agg ) { |agg, v| uv_info = edges[ [u,v] ]; agg += uv_info.last if (!uv_info.nil? && pi.pref(u,v) != uv_info.first); agg } }      
  end

  def sparse_local_improve(perm, edges)
    pi = perm
    is_moved = true
    while is_moved
      is_moved = false
      cost_pi = sample_weight_kt(pi, edges)
      (0..pi.size-1).each do |u|
        (0..pi.size-1).each do |v|
          npi = pi.move(u, v)
          cost_npi = sample_weight_kt(npi, edges)
          if cost_npi < cost_pi
            pi = npi 
            cost_pi = cost_npi
            is_moved = true           
          end
        end
      end    
    end
    pi
  end


  def vect_diff_zz(u,v)
    (u.keys - v.keys + v.keys).sort.inject({}) { |uv, fi| uv[fi.to_i] = (u[fi].nil? ? 0.0 : u[fi]) - (v[fi].nil? ? 0.0 : v[fi]); uv }
  end

  def vect_diff_fast_and_ugly2(u,v)
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

  def vect_diff_fast_and_ugly(u,v)
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



end 

#--
def self.constant_err_perm_whp( alternatives, oracle )
  Perm.new( (0..alternatives.size-1).to_a.shuffle.sort { |u,v| oracle.query( u, v ) == 1 ? -1 : 1 } )
end


end #module


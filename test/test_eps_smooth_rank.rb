require "minitest/unit"

require './'+File.dirname(__FILE__) + '/../lib/eps_smooth_rank.rb'
require './'+File.dirname(__FILE__) + '/../experiment/data_gen.rb'

class TestPerm < MiniTest::Unit::TestCase
    def setup
      @testArr = (0..99).to_a.shuffle
      @pi = EpsSmooth::Perm.new @testArr
    end

    def test_rank
      inx = rand 99
      assert_equal @pi.rank( inx ), @testArr[inx]
    end

    def test_inv_rank
      inx = rand 99
      assert_equal @pi.inv_rank( @testArr[inx] ), inx
    end

    def test_rank_inv_rank
      inx = rand 99
      assert_equal @pi.inv_rank( @pi.rank inx ), inx
    end

    def test_pref
      assert_equal @pi.pref( 0, 50 ), @testArr[0] < @testArr[50] ? 1 : 0
    end

    def test_size
      assert_equal @pi.size, @testArr.size
    end

    def test_move_first_to_last
      nPi = @pi.move @pi.inv_rank( 0 ), @pi.size-1
      (0..@pi.size-1).each do |u|
        assert_equal( nPi.rank( u ), @pi.rank( u )-1, "in between element failure" ) if @pi.rank( u ) > 0
        assert_equal( nPi.rank( u ), @pi.size-1, "last element failure" ) if @pi.rank( u ) == 0        
      end      
    end
    
    def test_move_high_inx_to_lower
      highInx  = @pi.size / 2
      lowInx = highInx / 2
      nPi = @pi.move @pi.inv_rank( highInx ), lowInx
      assert_equal nPi.rank( @pi.inv_rank( highInx ) ), lowInx
      (0..@pi.size-1).each do |u|
        assert_equal( nPi.rank( u ), @pi.rank( u )+1, "moving elements" ) if (lowInx..highInx-1).to_a.include? @pi.rank( u )
        assert_equal( nPi.rank( u ), @pi.rank( u ), "others" ) unless (lowInx..highInx).to_a.include? @pi.rank( u )
      end
    end

end



class TestSampler < MiniTest::Unit::TestCase
    def setup
      @data = gen_dist_prop_noise_data(100, 200, 0.1, 2)
    end

    def test_native_distributions
      pi = EpsSmooth::Perm.new [0,2,1,3]
      dists = EpsSrraNative::rank_pi_distributions( pi )
      assert dists.size == pi.size, "fail with dists size #{dists.size} / #{pi.size}"
      assert dists.first.size == pi.size, "fail with dists.first size #{dists.first} / #{pi}"
      dists.each do |d|
        assert_equal (0..d.size-1).inject(0.0){|s,i| s += d[i] - (i==0? 0.0 : d[i-1]) },1.0,"dist does not sum to 1.0 -- #{d}"
      end

      d = dists.last
      norm = [1.0/3, 1.0, 1.0/2, 0.0].inject(:+)
      d_should_be = [1.0/3/norm, 1.0/3/norm + 1.0/norm, 1.0/norm *( 1.0/3 + 1.0/2 + 1.0), 1.0/norm *(1.0/3 + 1.0/2 + 1.0)]      
      d.each_with_index { |v,i| assert(  (v-d_should_be[i]) < 0.0000001,"failed d content: #{d} / should be #{d_should_be}" ) }

      u_sample = EpsSrraNative::rank_sample_for_u( 0, dists.first, pi, 100)
    end

    def test_native_u_sample
      pi = EpsSmooth::Perm.new [0,2,1,3]
      dists = EpsSrraNative::rank_pi_distributions( pi )      

      u_sample = EpsSrraNative::rank_sample_for_u( 1, dists[1], pi, 10000)       
      u_sample.inject([0.0, 0.0, 0.0, 0.0 ]){ |m,r| m[ r[1] ] += 1.0/u_sample.size; m }.each_with_index do |epr,u|
        assert( ( epr - dists[1][u] + (u==0 ? 0 : dists[1][u-1]) ) < 0.02, "failed with conetent: u=#{u} -- estimation = #{epr} // real = #{dists[1][u]-(u==0 ? 0 : dists[1][u-1])} | all = #{dists[1]}" )
      end      
    end

    def test_native_pref_oracle
      oracle = EpsSmooth::PreferenceOracle.new @data["W"]
      pi = EpsSmooth::Perm.new (0..oracle.n-1).to_a.shuffle
      assert_equal oracle.kendalTauish_cost( pi ), oracle.native_kendalTauish_cost( pi )
    end

end

class Test_SVMLearner < MiniTest::Unit::TestCase
    def setup      
      @xs = (0..10).inject([]){|a,i| a << { 1=>rand, 2 => rand, 3=>rand }; a }
    end
    

    def test_svmrank
      ll = EpsSmooth::SVMLearner.new nil
      inxs = (0..@xs.size-1).to_a
      pairs = (0..@xs.size-1).to_a.map{ inxs.sample(2) }
      est = pairs.inject( { :labels =>[], :instances =>[], :weights =>[]} ) do |data,x|
        u1 = x[0]
        v1 =  x[1]      
        u2 = x[1]
        v2 =  x[0]
        pref_uv = rand > 0.5 ? 1 : 0
        
        data[:instances] << ll.vect_diff_fast_and_ugly( @xs[ u1 ], @xs[ v1 ] )
        data[:labels] << pref_uv
        data[:weights] << 1.0 / @xs.size

        data[:instances] << ll.vect_diff_fast_and_ugly( @xs[ u2 ], @xs[ v2 ] )
        data[:labels] << 1 - pref_uv
        data[:weights] << 1.0 / @xs.size

        data
      end

    end

end


MiniTest::Unit.autorun


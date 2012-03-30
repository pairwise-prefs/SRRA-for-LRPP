#include <ruby.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

static ID id_m_esrra;

static VALUE native_rank_pi_dist(VALUE m, VALUE pi){
  int size,u,v,i;
  float pr_uv;
  float **n_dist_all = NULL;
  float *n_norm = NULL;
  VALUE all_dist, dist_u;

  size = NUM2INT( rb_funcall(pi,rb_intern("size"), 0) );  

  /* mallocs */
  n_dist_all = malloc( size * sizeof(float*) );
  n_dist_all[0] = malloc( size*size*sizeof(float) ); /* multiple calls to malloc slows execution */
  for(i = 1; i < size; i++){
    n_dist_all[i] = n_dist_all[0] + i * size;
  }

  n_norm = (float*) malloc( size* sizeof(float) );
  /*..*/

  for(u=0; u<size; ++u){
    n_norm[u] = 0;
    for(v=0; v<size; ++v){
      pr_uv = ( v == u )? 0.f : 1.f / (float)fabs( FIX2INT(rb_funcall(pi,rb_intern("inv_rank"), 1, INT2FIX(u))) - FIX2INT(rb_funcall(pi,rb_intern("inv_rank"), 1, INT2FIX(v))) );
      n_dist_all[u][v] = (v == 0) ? 0.f : n_dist_all[u][v-1];
      n_dist_all[u][v] += pr_uv;
      n_norm[u] += pr_uv;
    }
  }

  all_dist = rb_ary_new();
  for(u=0; u<size; ++u){
    dist_u = rb_ary_new();   
    for(v=0; v<size; ++v){
      rb_ary_push( dist_u, rb_float_new( n_dist_all[u][v] / n_norm[u]) );
    }
    rb_ary_push( all_dist, dist_u );
  }

  /* free mem */
  free( (void*)n_dist_all[0] ); 
  free( (void*)n_dist_all );
  free( (void*)n_norm );
  /*..*/

  return all_dist;
}


/*...*/
static VALUE native_rank_pi_dist2(VALUE m, VALUE pi, VALUE beta){
  int size,u,v,i;
  float pr_uv;
  float **n_dist_all = NULL;
  float *n_norm = NULL;
  VALUE all_dist, dist_u;

  size = NUM2INT( rb_funcall(pi,rb_intern("size"), 0) );  

  /* mallocs */
  n_dist_all = malloc( size * sizeof(float*) );
  n_dist_all[0] = calloc( size*size, sizeof(float) ); /* multiple calls to malloc slows execution */
  for(i = 1; i < size; i++){
    n_dist_all[i] = n_dist_all[0] + i * size;
  }

  n_norm = (float*) malloc( size* sizeof(float) );
  /*..*/

  for(u=0; u<size; ++u){
    n_norm[u] = 0;
    for(v=0; v<size; ++v){
      pr_uv = ( v == u )? 0.f : 1.f / (float)pow( (float)fabs( FIX2INT(rb_funcall(pi,rb_intern("inv_rank"), 1, INT2FIX(u))) - FIX2INT(rb_funcall(pi,rb_intern("inv_rank"), 1, INT2FIX(v)))), NUM2DBL(beta));
      //      printf("--debug-- (u,v)=(%d,%d) pr_uv %.4f \n",u,v,pr_uv);
      n_dist_all[u][v] = (v == 0) ? 0.f : n_dist_all[u][v-1];
      n_dist_all[u][v] += pr_uv;
      n_norm[u] += pr_uv;
    }
  }

  all_dist = rb_ary_new();
  for(u=0; u<size; ++u){
    dist_u = rb_ary_new();   
    for(v=0; v<size; ++v){
      rb_ary_push( dist_u, rb_float_new( n_dist_all[u][v] / n_norm[u]) );
    }
    rb_ary_push( all_dist, dist_u );
  }

  /* free mem */
  free( (void*)n_dist_all[0] ); 
  free( (void*)n_dist_all );
  free( (void*)n_norm );
  /*..*/

  return all_dist;
}



/*...*/








/* generates an edge sample for alternative "u" (may include repititions)
   "m" is module name; "u" is the alternative id; "dist_arr" is a distribution array; "pi" is EpsSmooth::Perm instance; and "quota" defines the amount of edges to draw (may include repititions)
   OUTPUT: Array of Arrays. each element "x" has the format: "x[0]" - u; "x[1]" - v; "x[2]" - pi.pref(u,v); "x[3]" - sample-pr(u,v).
 */
static VALUE native_sample_for_u(VALUE m, VALUE u, VALUE dist_arr, VALUE pi, VALUE quota){
  float inv_max_rand = 1.f / (float) RAND_MAX;
  int i,int_v;
  float rand_pr;
  double pr_uv, weight_uv;
  VALUE v;
  VALUE all_sample = rb_ary_new();
  VALUE sample_info;

  static int do_srand = 1;
  if( do_srand ){
    //printf( "\n\n\n doing srand \n\n\n" );
    do_srand = 0;
    srand( (float)time(NULL) );
  }
  
  for(i=0; i<FIX2INT(quota); ++i){
    rand_pr = rand() * inv_max_rand;
    v = rb_funcall(dist_arr,rb_intern("dist_index"), 1, rb_float_new(rand_pr) );  // generalized binary search the index of rand_pr in "u's distribution"
    int_v = FIX2INT( v );

    //printf("\n\n rand_pr %f -- v: %d \n\n", rand_pr, int_v);

    sample_info = rb_ary_new(); // hold the sample info: u,v, pi.pref(u,v), sample-probability(u,v)
    rb_ary_push( sample_info, u ); // u
    rb_ary_push( sample_info, v ); // v
    rb_ary_push( sample_info, rb_funcall(pi,rb_intern("pref"),2, u,v ) ) ;  // pi.pref( u, v )

    pr_uv = ( NUM2DBL(rb_ary_entry(dist_arr, int_v )) - (int_v > 0 ? NUM2DBL(rb_ary_entry(dist_arr, int_v -1)) : 0.0 ) );
    rb_ary_push( sample_info, rb_float_new( pr_uv ) ); // sample pr(u,v)

    weight_uv = 1.0 / (NUM2DBL( quota ) * pr_uv); // (1.0 / |sample_size|*pr(u,v))
    rb_ary_push( sample_info, rb_float_new( weight_uv ) );
    
    rb_ary_push( all_sample, sample_info ); // aggregate info
  }
  return all_sample;
}


static VALUE native_aggregate_over_upperW(VALUE m, VALUE n){
  int u,v,int_n;
  double agg = 0.0;
  VALUE yield_result;

  if ( !rb_block_given_p() )
    rb_raise( rb_eArgError, "Expected a block");
  
  int_n = FIX2INT( n );

  for( u=0; u<int_n; ++u ){
    for( v=u+1; v<int_n; ++v ){
      yield_result = rb_yield_values( 2, INT2FIX(u), INT2FIX(v) );
      agg += NUM2DBL( yield_result );
    }
  }
    
  return rb_float_new( agg );
}


void esrra_native_library_initialize(){
  // rb_define_module_function(module, "function_name", implementation, number_of_args);
  rb_define_module_function(id_m_esrra, "rank_pi_distributions", native_rank_pi_dist, 1); 
  rb_define_module_function(id_m_esrra, "rank_pi_beta_distributions", native_rank_pi_dist2, 2); 
  rb_define_module_function(id_m_esrra, "rank_sample_for_u", native_sample_for_u, 4); 
  rb_define_module_function(id_m_esrra, "rank_aggergate_upperW", native_aggregate_over_upperW, 1);
}

/* initialize method for this extension. a ruby extension api requirement.
 */
void Init_eps_srra_native(){
  id_m_esrra = rb_define_module("EpsSrraNative");
  esrra_native_library_initialize();
}






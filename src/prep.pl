#!/usr/bin/perl
#__INSERT_LICENSE__

sub template_subst {
    my $text = shift();
    my %arg_list = %{shift()};
    my @new=();
    while ($text =~ /__(.*?)__/) {
	my $ident = $1;
	my $subst = $arg_list{$ident};
	if (! $subst) {
	    $subst = $ {$ident};
	}
	push @new, $`,$subst;
	$text = $';
    }
    $text= join("",@new,$text);
    return $text;
}

sub count_oper {
    my $ret="";
    while (my $oper = shift()) {
	$ret .= "op_count.$oper += cache->nelems;\n";
    }
    return $ret;
}

sub genone {
    $count_oper_size = 'cache->nelems';
    my %args = ('NAME' => shift(), 
		'OTHER_ARGS' => shift(), 
		'ELEM_OPERATIONS' => shift(),
		'COUNT_OPER' => copg(@_),
		'CACHE_OPERATIONS' => $cache_op);
    print template_subst($genone,\%args);
}

sub genone_all {
    genone('set','','**pto++ = **pfrom++') ;
    genone('add','','**pto++ += **pfrom++','sum') ;
    genone('rest','','**pto++ -= **pfrom++','sum') ;
    genone('mult','','**pto++ *= **pfrom++','mult') ;
    genone('div','','**pto++ /= **pfrom++','div') ;
    genone('rcp',',double c=1.',
	   '**pto++ = c/(**pfrom++)','div') ;
    genone('axpy',',double alpha',
	   '**pto++ += alpha * **pfrom++','mult','sum');
}


sub gen_setel {
    print template_subst($gensetel,{'NAME' => shift(), 
			    'OTHER_ARGS' => shift(), 
			    'ELEM_OPERATIONS' => shift(),
			    'CACHE_OPERATIONS' => $cache_op});
}

sub copg {
    my $ret="";
    while (my $oper = shift()) {
	$ret .= "op_count.$oper += $count_oper_size;\n";
    }
    return $ret;
}

sub gen_setel_all {
    $count_oper_size = '1';
    gen_setel('setel','','*cache->to = val');
    $COUNT_OPER = copg('sum');
    gen_setel('addel','','*cache->to += val','sum');
    $COUNT_OPER = copg('mult');
    gen_setel('multel','','*cache->to *= val','mult');
}

sub in_place {
    my %args = ('NAME' => shift(), 
		'ELEM_OPERATIONS' => shift(),
		'CACHE_OPERATIONS' => $cache_op,
		'FUN_ARGS' => 'const double val');
    my $other_args = shift();
    for my $oopt (keys %$other_args) { $args{$oopt} = $other_args->{$oopt}; }
    print template_subst($in_place,\%args);
}

sub in_place_all {
    in_place('set','**to++ = val',{'FUN_ARGS' => 'const double val=0.'});
    in_place('scale','**to++ *= val');
    in_place('add','**to++ += val');
    in_place('rcp','**to = val/(**to++)',{'FUN_ARGS' => 'const double val=1.'});
    in_place('fun','**to = (*fun_)(**to); **to++',
	     {'FUN_ARGS' => 'scalar_fun_t *fun_'});
    in_place('fun','**to = (*fun_)(**to,user_args); **to++',
  	     {'FUN_ARGS' => 'scalar_fun_with_args_t *fun_,void *user_args'});
}

sub gen_sum {
    my $args = {
	INI_LOOP => 'val=0', 
	NAME => shift(), 
	ELEM_OPERATIONS => shift(),
	COUNT_OPER => shift(),
	CACHE_OPERATIONS => $cache_op,
	OTHER_ARGS => '',
	OTHER_ARGS_P => '',
	POST_LOOP_OPS => '',
	PRE_LOOP_OPS => '',
	AFTER_ALL_STRIDES => '',
    };
    my $new_args = shift();
#    print "new gen_sum call:\n";
    while (my ($k,$v) = each %{$new_args}) {
#	print "arg: '$k' -> '$v'\n";
	$args->{$k} = $v;
    }
    # Put a comma for arguments
    $args->{'C'} = ($args->{'OTHER_ARGS'} !~ /^\s*$/ ? "," : "");
    print template_subst($gen_sum,$args);
}

sub gen_max {
    my $args = {
	'NAME' => shift(), 
	'INI_LOOP' => shift(), 
	'ELEM_OPERATIONS' => shift(),
	'COUNT_OPER' => shift(),
	'CACHE_OPERATIONS' => $cache_op,
	'PRE_LOOP_OPS'=>'double aux'
	};
    my $new_args = shift();
    while (my ($k,$v) = each %$new_args) { $args->{$k} = $v; }
    print template_subst($gen_sum,$args);
}

sub gen_sum_all {
    $count_oper_size = 'ntot';
    gen_sum('sum','val += **pa++',copg('sum'));
    gen_sum('sum_square','aux= **pa++; val += aux*aux',copg(qw(sum mult)),
	    { 'PRE_LOOP_OPS'=>'double aux'});
    gen_sum('sum_abs','val += fabs(**pa++)',copg(qw(sum abs)));
    gen_sum('norm_p','val += pow(fabs(**pa++),p)',copg(qw(sum abs)),
	    {'OTHER_ARGS'=>'const double p',
	     'OTHER_ARGS_P'=>'p',
	     'POST_LOOP_OPS'=>'val = pow(val,1./p)'});
    gen_sum('norm_p','val += int_pow(fabs(**pa++),p)',copg(qw(sum abs)),
	    {'OTHER_ARGS'=>'const int p',
	     'OTHER_ARGS_P'=>'p',
	     'POST_LOOP_OPS'=>'val = pow(val,1./double(p))'});
    gen_sum('assoc','f.set(f.fun2(**pa++,f.v()))','',
	    {PRE_LOOP_OPS=>'f.pre_all()',
             INI_LOOP => 'f.pre()', 
             POST_LOOP_OPS=>'f.post(); val = f.v()',
             AFTER_ALL_STRIDES=>'f.post_all()',
	     OTHER_ARGS=>'Fun2 &f',
	     OTHER_ARGS_P=>'f'});
    gen_max('max','val=**pa++;','aux=**pa++; if (aux>val) val=aux',copg(qw(sum fun)));
    gen_max('min','val=**pa++;','aux=**pa++; if (aux<val) val=aux',copg(qw(sum fun)));
    gen_max('max_abs','val=fabs(**pa++);',
	    'aux=fabs(**pa++); if (aux>val) val=aux',copg(qw(sum fun abs)));
    gen_max('min_abs','val=fabs(**pa++);',
	    'aux=fabs(**pa++); if (aux<val) val=aux',copg(qw(sum fun abs)));
}

sub export_vals_array{
    ## export to Newmat as right-value
    print template_subst($export_array,
			 {'CONST'=>'const',
			  'ARG'=>'Matrix & A',
			  'DEFINE_TARGET_POINTER'=>'double *to = A.Store()',
			  'DEFINE_NEWMAT_MAYBE'=>$DEFINE_NEWMAT}); 
    ## export to Newmat as left-value
    print template_subst($export_array,
			 {'CONST'=>'',
			  'ARG'=>'Matrix & A',
			  'DEFINE_TARGET_POINTER'=>'double *to = A.Store()',
			  'DEFINE_NEWMAT_MAYBE'=>$DEFINE_NEWMAT}); 
    ## export array as right-value
    print template_subst($export_array,
			 {'CONST'=>'const',
			  'ARG'=>'double *a',
			  'DEFINE_TARGET_POINTER'=>'double *to = a',
			  'DEFINE_NEWMAT_MAYBE'=>''}); 

    ## export array as left-value
    print template_subst($export_array,
			 {'CONST'=>'',
			  'ARG'=>'double *a',
			  'DEFINE_TARGET_POINTER'=>'double *to = a',
			  'DEFINE_NEWMAT_MAYBE'=>''}); 
}

/'/;
$warn_dont_modify =  <<EOT;
// DON'T EDIT MANUALLY THIS FILE !!!!!!
// This files automatically generated by ePerl from 
// the corresponding file (omitting the _eperl) name.
EOT
/'/;

1;

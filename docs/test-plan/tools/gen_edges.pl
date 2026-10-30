#!/usr/bin/perl
# Harper 2026-10-30: precision-matrix edge list v4 (Spence pmrev3 items 1-7 on top of v3).
use strict; use warnings;
my $T=$ENV{PM_TREE}//"$ENV{HOME}/harper-work/zf-precision-matrix";   # worktree; sources pinned to $SRC_COMMIT
my $IN="$T/docs/test-plan/inputs";
my $SRC_COMMIT="7955836a845bd6419dc41d21eda6e739e3b622c5";
my %PIN=("zf-inventory.csv"=>"ba79198eb3486944e3b43a898389b81a7ce8fb69e4f67450e49841e8fdac1d59","zb-inventory.csv"=>"b65c2bad80431500ab6280eb6ca545194407a35e9ae2f9719b5a6e5a38bdcd65","r4-required-results.csv"=>"2746dc4d05cca03d07db6cd8db600e86369a7f33b450c6b679b1d0c1bf0e8c6d","spence-gen_r4.py"=>"aa4d614617d086a49a2d625adbef1b516a7cac7ea37de58e391d40df2b51243f","FIXLIST.csv"=>"69391548209c157a3e8cc618abc3b474633b4d9cac89356f88ae28eb30d90887","fixes_after_1.0.12.txt"=>"03085a0acb3d5d0e51d8e8df13b3c4ad30d7843152d3d32024d6899c9f530c82","zf-audit-SUMMARY.md"=>"81cf6bb72636dfab03b2a25121724db4b37e3c16e0b896c30269db89be4122f1");
for my $f (sort keys %PIN){ my $o=`sha256sum "$IN/$f"`; my ($h)=$o=~/^([0-9a-f]{64})/; die "PIN MISMATCH $f: got ".($h//"missing")." want $PIN{$f}\n" unless ($h//"") eq $PIN{$f} }
system("git","-C",$T,"diff","--quiet",$SRC_COMMIT,"--","zenith-float-num/src","zenith-float-num/tests","zenith-bound/src","zenith-bound/tests")==0 or die "PIN MISMATCH: sources differ from $SRC_COMMIT\n";
my $INVZF="$IN/zf-inventory.csv"; my $INVZB="$IN/zb-inventory.csv"; my $R4="$IN/r4-required-results.csv";
my $OUT=shift or die "out"; my $N5OUT=shift or die "n5 out"; my $SHOUT=shift or die "shards out";

sub parse_csv { my ($txt)=@_; my @rows; my @f; my $fld=''; my $q=0; my $i=0; my $n=length $txt;
  while($i<$n){ my $c=substr($txt,$i,1);
    if($q){ if($c eq '"'){ if(substr($txt,$i+1,1) eq '"'){$fld.='"';$i++} else {$q=0} } else {$fld.=$c} }
    else { if($c eq '"'){$q=1} elsif($c eq ','){push @f,$fld;$fld=''} elsif($c eq "\n"){push @f,$fld;push @rows,[@f];@f=();$fld=''} elsif($c ne "\r"){$fld.=$c} }
    $i++ } push(@f,$fld),push(@rows,[@f]) if $fld ne '' || @f; return \@rows }
sub slurp { local $/; open my $h,'<',$_[0] or die $_[0]; my $s=<$h>; $s }
sub sha { my $o=`sha256sum "$T/$_[0]"`; ($o=~/^([0-9a-f]{64})/)?$1:'MISSING' }
# ---------- families (R1: Type::name)
my @fam;  # [crate,key,bare,path,sig,category,acc]
{ my $r=parse_csv(slurp($INVZF)); my @h=@{shift @$r}; my %seen;
  for my $x (@$r){ my %r; @r{@h}=@$x; next unless ($r{reachable_from_crate_root}//'') =~ /^(yes|macro-impl)$/; next if ($r{status}//'')=~/trait-decl/;
    my $ic=$r{impl_context}//''; $ic=~s/^\s+|\s+$//g; my $f=$r{path}; $f=~s/:.*//; $f=~s{^zenith-float-num/src/}{};
    my $ty; if($ic=~/^impl(?:<[^>]*>)?\s+(?:[\w:<>, ]+\s+for\s+)?(\$?[A-Za-z_]\w*)/){ $ty=$1; $ty="$ty\[$f]" if $ty=~/^\$/ } else { ($ty=$f)=~s/\.rs$//; $ty=~s{/}{::}g }
    my $k="$ty\::$r{name}"; next if $seen{$k}++; my $acc=$r{documented_accuracy}//''; $acc=~s/\s+$//;
    my $cat=$r{category}; $cat='random-feature' if $r{path}=~m{common/test_rng\.rs};
    push @fam,['zf',$k,$r{name},$r{path},$r{signature},$cat,$acc] } }
push @fam,['zf','common::test_rng::DEFAULT_RANDOM_SEED','DEFAULT_RANDOM_SEED','zenith-float-num/src/lib.rs:204-207','pub const DEFAULT_RANDOM_SEED: u64','random-feature',''];
{ my $txt=slurp($INVZB); $txt=(split /# SECTION 2/,$txt)[0]; $txt=join "\n", grep {!/^#/} split /\n/,$txt; $txt.="\n";
  my $r=parse_csv($txt); my @h=@{shift @$r}; my %seen;
  for my $x (@$r){ my %r; @r{@h}=@$x; next unless $r{name}; my $k="$r{type}::$r{name}"; next if $seen{$k}++;
    push @fam,['zb',$k,$r{name},"zenith-bound/$r{path}",$r{signature},"zb-$r{type}",''] } }
# ---------- R4 (Spence) rows for zf category=other, keyed by bare name
my %r4; { my $r=parse_csv(slurp($R4)); shift @$r; for(@$r){ $r4{"$$_[0]|$$_[1]"}=[$$_[2],$$_[3],$$_[4]] } }
my %r4fam=map {(split /\|/)[0]=>1} keys %r4;
# item 1: no "or documented" anywhere
for my $k (keys %r4){ for (@{$r4{$k}}){ s/ or documented NaN//g; s/ or documented range//g } }
# Spence G1-G3 (spence-gen_r4.py lines 6-19, copied to inputs/)
my %G123=map {$_=>1} qw(re im mid rad num den coeff coeffs degree shape radix precision emax emin err inexact inexact_accumulated classify as_raw_parts significant_bits cached_bit_len cache_info inner value source deref deref_mut clone default nan zero one i with with_radix try_new filled filled_2d reshape transpose conj inv_sign soft_neg map_exact invert gcd mod_pow mod_inv miller_rabin pollard_rho soft_cmp compare contains commensurable_shift);
# item 5: R4 texts restored for these names (all value classes)
my %R4X=(as_raw_parts=>'exact and mode-independent: returns the stored (sign, exponent, limbs, precision) unchanged; identical for every rm and on both targets',
 significant_bits=>'exact and mode-independent: the exact count of significant bits of the stored mantissa; identical for every rm and on both targets',
 random_seed=>'determined by the seed: output is a pure function of (seed, p); digest-pinned in the gold; byte-identical on i686 and x86_64',
 reseed_random=>'determined by the seed: after reseed_random(s) the stream equals seeded_random(s); digest-pinned in the gold; byte-identical on i686 and x86_64',
 seeded_random=>'determined by the seed: output is a pure function of (seed, p); digest-pinned in the gold; byte-identical on i686 and x86_64');
# item 12: Rule R pins (tests/nan_error_kinds.rs EXPECTED)
my %NANPIN; { open my $nf,'<',"$T/zenith-float-num/tests/nan_error_kinds.rs" or die; my $ln=0; while(<$nf>){ $ln++; if(/^(\S+) (AB|BA|AN|NA) => (.*)$/){ $NANPIN{"$1|$2"}="$3 (pinned: zenith-float-num/tests/nan_error_kinds.rs:$ln)" } } }
my %SZ=(
 si=>{r=>'si(+0)=+0, si(-0)=-0 (odd)', c=>'si(+-0+-0i)=+-0+-0i with both signs kept (odd, conj-symmetric); si(-1+0i)=-Si(1)+0i, si(-1-0i)=-Si(1)-0i'},
 ci=>{r=>'ci(+0)=-Inf and ci(-0)=-Inf, each with error kind Error::DivisionByZero', c=>'ci(+0+0i)=-Inf+0i, ci(+0-0i)=-Inf-0i, ci(-0+0i)=-Inf+i*pi, ci(-0-0i)=-Inf-i*pi (each with Error::DivisionByZero); cut: ci(-1+0i)=Ci(1)+i*pi, ci(-1-0i)=Ci(1)-i*pi (DLMF 6.4.4)'},
 li=>{r=>'li(+0)=+0 and li(-0)=+0 exactly, no error kind (Spence-approved li(+0)=+0)', c=>'li(+0+0i)=+0+0i, li(+0-0i)=+0-0i, li(-0+0i)=+0+0i, li(-0-0i)=+0-0i (no error kind); cut x<0: li(-1+0i)=conj(li(-1-0i)), imaginary part +pi side for +0i and -pi side for -0i'},
 nth_root=>{r=>'nth_root(+0,n)=+0; nth_root(-0,n)=-0 for odd n and +0 for even n (n>0); n<0: +0 gives +Inf, -0 gives -Inf (odd n) or +Inf (even n), with Error::DivisionByZero (IEEE 754-2019 rootn)', c=>'nth_root(+-0+-0i,n)=+0+-0i (sign of im kept, as csqrt C11 G.6.4.2); cut: nth_root(-1+0i,n)=cos(pi/n)+i*sin(pi/n), nth_root(-1-0i,n)=cos(pi/n)-i*sin(pi/n)'},
 arg=>{r=>'arg(+0)=+0, arg(-0)=+pi', c=>'arg = atan2(im,re) (C11 7.3.9.1 carg, IEEE 754-2019 9.2.1): arg(+0+0i)=+0, arg(+0-0i)=-0, arg(-0+0i)=+pi, arg(-0-0i)=-pi, arg(-1+0i)=+pi, arg(-1-0i)=-pi'});
my %KN=(A=>'Error::InvalidArgument',B=>'Error::DivisionByZero',N=>'no kind (None)');
# item 3: N=5 spot-check list (wrong answer ever), citations into docs/test-plan/inputs/
my $FX='docs/test-plan/inputs/FIXLIST.csv'; my $SU='docs/test-plan/inputs/zf-audit-SUMMARY.md'; my $FA='docs/test-plan/inputs/fixes_after_1.0.12.txt';
my @N5=( # [crate, type-regex, bare-names, citation]
 ['zf','.*',[qw(sin)],"$FX:15 (ZF-01); $FX:127 (ZF-M22); $FA:87"],
 ['zf','.*',[qw(cos sin_cos tan)],"$FX:127 (ZF-M22 all trig_reduce_extra users); $FA:87 (cos, sin_cos); $FA:98 (trig)"],
 ['zf','.*',[qw(rem_pi)],"$FX:17 (ZF-03 rem_pi rounding); $FX:31 (ZF-17); $FX:127"],
 ['zf','^ExactComplex$',['*'],"$FX:16 (ZF-02 complex word-rounded); $SU:29 (#8 complex mul/div); $FA:42-44"],
 ['zf','.*',[qw(gamma)],"$FX:17 (ZF-03); $SU:24 (#3); $FA:31; $FA:35; $FA:98"],
 ['zf','.*',[qw(ln_gamma)],"$SU:24 (#3); $FA:38"],
 ['zf','.*',[qw(digamma)],"$SU:24 (#3); $FA:30; $FA:34"],
 ['zf','.*',[qw(exp)],"$FX:17 (ZF-03 exp far-underflow); $FA:87; $FA:98"],
 ['zf','.*',[qw(gamma_pdf beta_pdf poisson_pmf binomial_pmf chi_squared_cdf student_t_pdf)],"$FX:19 (ZF-05)"],
 ['zf','.*',[qw(add)],"$FX:20 (ZF-06); $FX:117 (ZF-M19 add panics on i686); $SU:25 (#4); $FA:101"],
 ['zf','.*',[qw(sub mul div)],"$FX:20 (ZF-06); $SU:25 (#4); $FA:101; $FA:44 (div)"],
 ['zf','.*',[qw(fma)],"$SU:25 (#4); $FA:45; $FA:83"],
 ['zf','.*',[qw(reciprocal)],"$SU:25 (#4); $FA:19"],
 ['zf','.*',[qw(pow powi powsi)],"$SU:22 (#1); $FA:8 (powsi); $FA:21-25; $FA:57-62"],
 ['zf','.*',[qw(round)],"$FA:21-25"],
 ['zf','.*',[qw(riemann_zeta)],"$SU:23 (#2); $FA:37; $FA:106"],
 ['zf','.*',[qw(expm1)],"$SU:26 (#5); $FA:33; $FA:87"],
 ['zf','.*',[qw(log1p)],"$SU:26 (#5); $FA:40; $FA:86; $FA:87"],
 ['zf','.*',[qw(log10)],"$SU:26 (#5); $FA:52"],
 ['zf','.*',[qw(log2)],"$SU:26 (#5); $FA:15"],
 ['zf','.*',[qw(log ln)],"$SU:26 (#5 log-with-base); $FA:18; $FA:98"],
 ['zf','.*',[qw(acosh)],"$SU:26 (#5); $FA:49"],
 ['zf','.*',[qw(atan2)],"$SU:26 (#5); $FA:46"],
 ['zf','.*',[qw(ldexp)],"$SU:26 (#5); $FA:48"],
 ['zf','.*',[qw(nth_root)],"$SU:26 (#5); $FA:50-51"],
 ['zf','.*',[qw(pi e ln_2 ln_10)],"$SU:27 (#6 constants); $FA:53"],
 ['zf','.*',[qw(erf erfc)],"$SU:28 (#7); $FA:39; $FA:41; $FA:98"],
 ['zf','.*',[qw(sinh cosh tanh)],"$FA:87; $FA:98 (hyperbolic)"],
 ['zf','.*',[qw(sqrt)],"$FA:88; $FA:99"],
 ['zf','.*',[qw(cbrt)],"$FA:98"],
 ['zf','.*',[qw(jacobi_am)],"$FA:89"],
 ['zf','.*',[qw(bessel_j0 bessel_j1)],"$FA:36"],
 ['zf','.*',[qw(ei ci li)],"$FA:98"],
 ['zf','.*',[qw(elliptic_k elliptic_e)],"$FA:98"],
 ['zf','.*',[qw(clausen_cl2 clausen_cl3 weber_e barnes_g)],"$FA:106"],
 (map { my ($n,$s,$d)=@$_; [q(zf),q(.*),[$n],"$FX:18 (ZF-04 S1 wrong-answer, 33 claimed); zenith-float-num/src/$s at SRC_COMMIT 7955836 (doc: Accuracy: unspecified)"] } ([qw(scorer_gi catalog.rs:575-576)],[qw(scorer_hi catalog.rs:603-604)],[qw(kelvin_ber_nu catalog.rs:629-630)],[qw(kelvin_bei_nu catalog.rs:658-659)],[qw(kelvin_ker_nu catalog.rs:687-688)],[qw(kelvin_kei_nu catalog.rs:722-723)],[qw(kelvin_ber catalog.rs:757-758)],[qw(kelvin_bei catalog.rs:772-773)],[qw(kelvin_ker catalog.rs:787-788)],[qw(kelvin_kei catalog.rs:802-803)],[qw(struve_h catalog.rs:821-822)],[qw(anger_j catalog.rs:916-917)],[qw(weber_e catalog.rs:962-963)],[qw(clausen_cl2 catalog.rs:1007-1008)],[qw(clausen_cl3 catalog.rs:1038-1039)],[qw(hurwitz_zeta catalog.rs:1066-1067)],[qw(riemann_zeta catalog.rs:1115-1116)],[qw(polygamma catalog.rs:1148-1149)],[qw(barnes_g catalog.rs:1181-1182)],[qw(hypergeom_0f1 catalog.rs:1328-1329)],[qw(hypergeom_1f1 catalog.rs:1347-1348)],[qw(hypergeom_pfq catalog.rs:1389-1390)],[qw(lambert_w0 catalog.rs:1469-1470)],[qw(lambert_wm1 catalog.rs:1484-1485)],[qw(polylog catalog.rs:1581-1582)],[qw(meijer_g heun_fox.rs:576-577)],[qw(fox_h heun_fox.rs:640-641)],[qw(heun_g heun_fox.rs:771-772)],[qw(heun_c heun_fox.rs:967-968)])), (map { [q(zf),q(.*),[$_->[0]],"$FX:18 (ZF-04 range catalog.rs:575-1581); zenith-float-num/src/$_->[1] at SRC_COMMIT 7955836 (doc: Accuracy: correctly-rounded, NOT unspecified; outside the 29 unspecified, kept at N=5 because it is in the cited range and goes through round_exact)"] } ([qw(jacobi_arcsn catalog.rs:1254-1255)],[qw(jacobi_arccn catalog.rs:1273-1274)],[qw(jacobi_arcdn catalog.rs:1300-1301)])),
 ['zb','^Ball$',[qw(add sub mul_exact powi sqrt div inv hypot abs)],"$FX:67 (ZB-01 eps radius 0); $FX:68 (ZB-02 same 9 ops)"],
);
my (%n5,%n5hit);

# ---------- golds already in tree (pins for rows whose bare name appears in them)
my %goldfor; for my $g (glob "$T/zenith-float-num/tests/oracle/goldens/*.jsonl"){ my $rel=$g; $rel=~s{^\Q$T\E/}{}; my $s=slurp($g);
  my %n; $n{lc $1}=1 while $s=~/"(?:op|fn|c)"\s*:\s*"([A-Za-z_0-9]+)"/g; my $h=sha($rel); push @{$goldfor{zf}{$_}},"$rel\@sha256:$h" for keys %n }
for my $g (glob("$T/zenith-bound/tests/goldens/*.txt"), glob("$T/zenith-bound/tests/goldens/slices/*/*.txt")){ my $rel=$g; $rel=~s{^\Q$T\E/}{}; my $h=sha($rel);
  open my $fh,'<',$g; my %n; while(<$fh>){ next if /^#/; my @w=split; next unless @w>1; $n{$w[0]}=1; $n{$w[1]}=1 if $rel=~/matrix_rows/ } push @{$goldfor{zb}{$_}},"$rel\@sha256:$h" for keys %n }
my $MPFR='MPFR 4.2.2 + GMP 6.3.0 via rug 1.20.1 / gmp-mpfr-sys 1.6.8 (Cargo.lock @7955836), feature mpfr-tests, live on x86_64-unknown-linux-gnu';
my $ENV='~/team-env/oracle-2026-10-30 via ORACLE_PY (mpmath 1.3.0, sympy 1.14.0, gmpy2 2.3.2, python-flint 0.9.0)';
my %mpfrnative=map {$_=>1} qw(add sub mul div sqrt cbrt rec_sqrt recip sqr exp exp2 exp10 expm1 ln log log2 log10 ln_1p log1p sin cos tan sec csc cot asin acos atan atan2 sinh cosh tanh sech csch coth asinh acosh atanh gamma ln_gamma lgamma digamma zeta erf erfc j0 j1 jn y0 y1 yn ai agm hypot pow powi powsi fma fms li2 eint fmod remainder rint round floor ceil trunc frac abs neg sin_cos sinh_cosh nth_root beta gamma_inc);
my %trig=map {$_=>1} qw(sin cos tan sec csc cot sin_cos sincos rem_pi);
my %i64x=map {$_=>1} qw(scalb ldexp ilogb frexp);
my $W=0;
my @MODES=qw(None Up Down ToZero FromZero ToEven ToOdd);
my %mpfrmode=(None=>'RNDZ at rp(p,64) (R6)',Up=>'RNDU',Down=>'RNDD',ToZero=>'RNDZ',FromZero=>'RNDA',ToEven=>'RNDN',ToOdd=>'RNDN, exact tie flipped to the odd neighbour by the oracle (tie detected as exact at p+1 bits)');
my @out; my %cnt; our @PROOFS;
sub row { my %a=@_; push @out,[@a{qw(crate family source category group class class_def inputs modes tier required_result refusal_variant oracle_x86_64 oracle_i686 check verdict_space proof)}]; $cnt{"$a{crate} $a{verdict_space}"}++ }
for my $F (@fam){ my ($crate,$key,$bare,$src,$sig,$cat,$acc)=@$F; $sig//='';
  my ($ty)=$key=~/^(.*)::[^:]+$/; $ty//='';
  for my $i (0..$#N5){ my ($c,$tr,$names,$cit)=@{$N5[$i]}; next unless $c eq $crate && $ty=~/$tr/; if(grep {$_ eq '*' || $_ eq $bare} @$names){ push @{$n5{$key}},$cit; $n5hit{$i}{$bare}=1 } }
  my $hasp = $sig=~/\b(p|prec|precision)\s*:\s*(usize|Precision)\b/ ? 1:0;
  my $hasrm= $sig=~/RoundingMode/ ? 1:0;
  my $ret = $sig=~/->\s*Option</ ? 'option' : $sig=~/->\s*(?:core::result::)?Result<|->\s*BoundResult</ ? 'result' : $sig=~/->\s*bool\b/ ? 'bool' : 'value';
  my ($args)=$sig=~/\((.*)\)/s; $args//=''; my $nops=()=$args=~/(?:&\s*(?:mut\s+)?self\b|\bself\b|:\s*&?\s*(?:Self|ExactNum|ExactComplex|ExactInt|Ball|CBall|Mag)\b)/g;
  my $binary = $nops>=2 ? 1:0;
  my $wp = $acc=~/working-precision, 32 bits/ ? 1:0;
  my $other = ($crate eq 'zf' && $cat eq 'other') ? 1:0;
  my $grp = $crate eq 'zb' ? 'zb-enclosure' : $cat eq 'random-feature' ? 'G13 random (feature=random, R2)' : $other ? 'R4' : $cat=~/^(arith|elementary|special|trig|constant|rounding)$/ ? 'numeric-CR' : $cat=~/^(conversion|io\/format)$/ ? 'conversion' : $cat eq 'compare' ? 'compare' : 'numeric-CR';
  my $V = sub { my $v=shift; return "Err(BoundError::$v)" if $crate eq 'zb'; return 'None' if $ret eq 'option'; return "Err(Error::$v)" if $ret eq 'result'; return "false (documented refusal; $v-class input)" if $ret eq 'bool'; "NaN whose error kind is Error::$v" };
  my $proofp="docs/test-plan/pm-proofs/$crate/".($key=~s/[^A-Za-z0-9_:.-]/_/gr).'.md';
  my $gold = sub { my $g=($goldfor{$crate}{lc $bare}||[]); @$g ? join(' | ',map {"gold:$_"} @$g) : "gold:".($crate eq 'zb'?'zenith-bound/tests/goldens/pm':'zenith-float-num/tests/pm/golds')."/$crate-".($key=~s/[^A-Za-z0-9_]+/_/gr).".jsonl\@sha256:UNGENERATED (pm-coverage FAIL until written by PM-3 and pinned in SHA256SUMS)" };
  my $ox = $crate eq 'zb' ? "goldlib (zenith-bound/tools, ORACLE_PY=$ENV; python-flint 0.9.0 arb at prec+256), ".$gold->()
        : ($mpfrnative{$bare} && $key=~/^ExactNum::/) ? $MPFR
        : ($grp eq 'G13 random (feature=random, R2)') ? 'digest:'.$gold->()
        : ($bare=~/gcd|lcm|mod_|egcd|pollard|miller|factor|binomial|isqrt|int/i) ? "python-flint 0.9.0 exact in $ENV, ".$gold->()
        : "mpmath 1.3.0 Ziv pair p+256/p+512 must agree, in $ENV, ".$gold->();
  my $oi = $crate eq 'zb' ? 'same pinned gold file + SHA as x86_64; byte-identical BallRow digest to the x86_64 run'
        : 'gold:zenith-float-num/tests/pm/golds/i686/'.($key=~s/[^A-Za-z0-9_]+/_/gr).".jsonl\@sha256:UNGENERATED (cut from the MPFR/gold-checked x86_64 run at PM-3; pm-coverage FAIL until pinned)";
  my $chk = $crate eq 'zb' ? 'containment: exact value in [mid-rad, mid+rad]; rad <= tightness baseline' : 'exact limbs + sign + inexact flag + error kind equal the pinned oracle';
  my %base=(crate=>$crate,family=>$key,source=>$src,category=>$cat,group=>$grp,oracle_x86_64=>$ox,oracle_i686=>$oi,check=>$chk,proof=>'',refusal_variant=>'');
  my $cr = $crate eq 'zb' ? 'enclosure containing the exact value (sound), rad <= baseline' : $wp ? 'R3: correctly rounded at the requested p in the requested rm, or refuse with '.$V->('InvalidArgument').'; a result accurate only to 32 bits is FAIL (FIXLIST: working-precision families)' : 'correctly rounded at p in rm';
  my $R = sub { my ($cls,$def,$inp,$modes,$tier,$req,$var,%o)=@_;
    if($crate eq 'zf' && $R4X{$bare} && !$o{keep} && !$o{na} && $cls!~/refuse|^nan/){ $req="R4 $R4X{$bare}" }
    elsif($other && !$o{keep} && !$o{na}){ my $pc=$o{r4class}//$cls; my $h=$r4{"$bare|$pc"}; if($h){ $req="R4 $$h[0]"; $var=$$h[1] if !$var; $req.=" [convention: $$h[2]]"; $req=~s/ or documented (NaN|range)//g } else { $req="R4 MISSING for $bare|$pc"; } }
    my $vs = $o{vs} // ($o{na} ? 'NA_PROOF' : 'PASS|FAIL'); $req=~s/ or documented (NaN|range)//g;
    if($o{g123}){ push @PROOFS,[$proofp,$key,$src,$bare] }
    row(%base,class=>$cls,class_def=>$def,inputs=>$inp,modes=>$modes,tier=>$tier,required_result=>$req,refusal_variant=>$var//'',verdict_space=>$vs,proof=>$o{na}?"$proofp#$cls (committed proof, signed by a second agent; pm-coverage FAIL if missing)":'', ($o{na}?(oracle_x86_64=>'n/a (proof row)',oracle_i686=>'n/a (proof row)',check=>'proof file present, cites signature/type at '.$src.', second-agent sign-off line present'):()) ) };
  my $SKVS='PASS|FAIL|SKIPPED_AFTER_TIMEOUT'; my $SKIP='MAX-timeout skip (Paul, Spence pmrev3 6): the x86_64 MAX row of a family runs first; if it times out, the same family\'s i686 MAX row is not run and is recorded SKIPPED_AFTER_TIMEOUT, which counts as FAIL in the verdict space and in pm-coverage and is never omitted; the family\'s RC stays blocked until its FIXLIST fix lands, and every skipped row then reruns in full';
  my $M7='ToEven at p in {32,64,128,256,512,1024}; heavy: all 7 modes x axis';
  my $X7='fast: all 7 modes at p=64; heavy: all 7 modes x axis';
  my $FX='fast+heavy'; my $NA=!$hasrm;
  my $rmmodes = $hasrm ? $X7 : 'n/a: no RoundingMode parameter';
  # value classes
  $R->('normal','random normal operands in domain','splitmix seeded, 16 per cell',$M7,$FX,$cr);
  $R->('negative','negative operands','splitmix seeded',$M7,$FX,$cr);
  if($crate eq 'zf' && $bare eq 'ci'){ $R->('zero','+0 operand','+0 (and 0+0i for ExactComplex)',$M7,$FX,'Ci(+0) = -Inf with error kind Error::DivisionByZero (log singularity, DLMF 6.6.2); the same in every rm and on both targets','Error::DivisionByZero (-Inf)',keep=>1) }
  elsif($crate eq 'zf' && $bare eq 'li'){ $R->('zero','+0 operand','+0 (and 0+0i for ExactComplex)',$M7,$FX,'li(+0) = +0 exactly (DLMF 6.2.8: li(x)=Ei(ln x), Ei(-Inf)=0); no error kind; the pole li(1) is the out-of-domain row','',keep=>1) }
  else { $R->('zero','+0 operand(s)','+0',$M7,$FX,$crate eq 'zb'?'enclosure; zero mid is +0 (signed zero normalised, PLAN L1)':$cr) }
  # item 6: operand classes restored
  $R->('huge-exp','operand(s) with exponent in [EXPONENT_MAX-64, EXPONENT_MAX] (huge operand, not a huge result)','constructed from exact exponents',$hasrm?$X7:$M7,$FX,$crate eq 'zb'?'enclosure, or the all-reals ball (rad=Infinity) when the midpoint overflows (CONTRACT 10.4)':'correct bits at p in rm, else Error::ExponentOverflow(sign) per mode (away from zero -> +-Inf, toward zero -> max finite)',$crate eq 'zb'?'':'ExponentOverflow(Sign)');
  $R->('tiny-exp','operand(s) with exponent in [EXPONENT_MIN, EXPONENT_MIN+64] (tiny operand)','constructed from exact exponents',$hasrm?$X7:$M7,$FX,$crate eq 'zb'?'enclosure; the radius absorbs any lost value and is never 0 (ZB-01, CONTRACT 3.2)':'correct bits incl. subnormal results at p in rm, else the signed zero of the exact sign per mode','');
  $R->('subnormal-like','unnormalised raw parts: from_raw_parts with leading zero limbs / top bit clear (R4 subnormal-like)','ExactNum::from_raw_parts(limbs with top bit 0 or leading zero words, exp, sign, p)',$hasrm?$X7:$M7,$FX,$crate eq 'zb'?'enclosure built from the normalised midpoint (bridge.rs); an unrepresentable midpoint is refused with Err(BoundError::Representation)':'from_raw_parts normalises exactly and the family result is then correct at p in rm, or from_raw_parts refuses with Error::InvalidArgument (from_raw_parts checks); never a value computed from the unnormalised mantissa',$crate eq 'zb'?'Err(BoundError::Representation)':'Error::InvalidArgument');
  if($crate eq 'zb' && $bare eq 'div'){
    $R->('divisor-contains-0','divisor ball with |m_b| <= r_b and r_b > 0','m_b in {0, ±r_b/2, ±r_b}, r_b > 0','n/a','fast+heavy','0 in the divisor with rad>0 gives all reals','Ok(all-reals ball: mid=+0, rad=Infinity)');
    $R->('divisor-exact-zero','divisor mid 0, rad 0','b = 0 ± 0','n/a','fast+heavy','an exact 0 divisor gives Domain',"Err(BoundError::Domain)");
  }
  if($crate eq 'zf' && $SZ{$bare}){ my $cx=$key=~/^ExactComplex::/?'c':'r'; my $t=$SZ{$bare}{$cx}; $R->('signed-zero',"signed zeros and cut points ".($cx eq 'c'?'(complex: every sign of re and im zero)':'(real)'),$cx eq 'c'?'+0+0i, +0-0i, -0+0i, -0-0i and the cut points listed':'+0, -0','fast: all 7 modes at p=64; heavy: all 7 modes x axis',$FX,"exact expected values (C11 Annex G conjugate symmetry f(conj z)=conj f(z), principal branch, the sign of a zero imaginary part chooses the side of the cut; e.g. ln(-1-0i)=-i*pi): $t; values are exact or correctly rounded at p in rm; any mismatch is a red row (FIXLIST doc row for LIBRARY.md:399 is Paul's)",($t=~/DivisionByZero/?'Error::DivisionByZero':''),keep=>1) } else {
  $R->('signed-zero','-0 and +0 operands','-0, +0 (and pairs for binary)',$hasrm?$X7:$M7,$FX,$crate eq 'zb'?'enclosure; -0 mid normalised to +0':'correct value and IEEE/MPFR sign of zero in every mode (exact cancel gives -0 only in Down)'); }
  $R->('inf','±Inf operand(s)','±Inf',$M7,$FX,$crate eq 'zb'?'refusal':'correct limit (IEEE 754-2019 7.2/9.2, MPFR semantics) or refusal',$crate eq 'zb'?$V->('Representation'):$V->('InvalidArgument').' where the limit is undefined');
  my $cplx = $key=~/^ExactComplex::/ ? 1:0;
  if($crate eq 'zf' && ($binary || $cplx)){ for my $pp (qw(AB BA AN NA)){ my ($k1,$k2)=split //,$pp;
      my $pinkey = $cplx ? ($bare=~/^(asinh|ln_gamma)$/ ? $bare : "c$bare") : $bare; my $pin=$NANPIN{"$pinkey|$pp"};
      my ($def,$req);
      if($cplx && !$binary){ $def="complex NaN operand z with kinds (re,im) = ($k1,$k2)"; $req= $pin ? "inherited: $pin" : $bare eq 'asinh' ? '' : "inherited: re=$KN{$k1}, im=$KN{$k1} (the re part's kind wins wherever the parts mix, as pinned for cexp/csqrt/asinh; Rule R, CHANGELOG.md:131-132)" }
      elsif($cplx){ $def="complex NaN operands a, b with all parts of kind $k1 (a) and $k2 (b)"; $req= $pin ? "inherited: $pin" : "inherited: re=$KN{$k1}, im=$KN{$k1} (componentwise Rule R: the first operand's kind wins per part; CHANGELOG.md:131-132)" }
      else { $def="NaN operands (first,second) kinds = ($k1,$k2); A=InvalidArgument B=DivisionByZero N=no kind"; $req= "inherited: $KN{$k1} from the first operand (Rule R: first operand's kind wins, including kind-less N)".($pin?"; $pin":"; Rule R stated at zenith-float-num/tests/nan_error_kinds.rs:4-7 and CHANGELOG.md:131-132") }
      $req="[PM1-UNPINNED] required result: FAIL until this case is pinned in zenith-float-num/tests/nan_error_kinds.rs EXPECTED at PM-1 (a run whose kind matches the expectation below still FAILs until then); expectation to pin: $req" if !$pin;
      $R->("nan-$pp",$def,'NaN operands with pinned kinds',$M7,$FX,$req,$k1 eq 'N'?'no kind (None)':$KN{$k1},keep=>1) } }
  elsif($crate eq 'zf'){ for my $k (qw(A B N)){ $R->("nan-$k","NaN operand with kind $k (A=InvalidArgument B=DivisionByZero N=no kind)",'NaN operand with pinned kind',$M7,$FX,"NaN result, inherited kind $KN{$k} (unary propagation keeps the operand's kind; nan_error_kinds.rs Rule R)",$k eq 'N'?'no kind (None)':$KN{$k},keep=>1) } }
  else { $R->('nan','NaN operand','NaN midpoint (kind A, B and N)',$M7,$FX,'refusal (NaN midpoint never stored): every operand kind maps to the same variant',$V->('Representation')) }
  if($crate eq 'zf' && $bare=~/^(ci|li)$/){ my $pt=$bare eq 'ci'?'Ci(+0) and Ci(-0) (real: -0 is on the cut)':'li(1)';
    $R->('out-of-domain',"pole $pt; cut points",'exact pole; negative reals on the cut with +-0 imaginary part',$M7,$FX,($bare eq 'ci'?'Ci(0)':'li(1)').' = -Inf with error kind Error::DivisionByZero (logarithmic pole, DLMF 6.2.8/6.6.2); identical in every rm and on both targets; cut points: principal value with +0 imaginary side','Error::DivisionByZero (-Inf)',keep=>1) }
  else {
  $R->('out-of-domain','argument outside the mathematical domain, at a pole, or on a cut','domain: per family table; pole: exact pole',$M7,$FX,$crate eq 'zb'?'Strict: Domain if any point is outside the domain; Clip per CONTRACT 10.2/10.2a':'refusal with the real variant: domain -> InvalidArgument, pole -> DivisionByZero, cap/budget -> PrecisionRetryExhausted',$crate eq 'zb'?$V->('Domain'):join('; ',$V->('InvalidArgument'),$V->('DivisionByZero'),$V->('PrecisionRetryExhausted'))); }
  $R->('cancellation','x - x, near-cancelling pairs','x and nextafter(x) pairs',$hasrm?$X7:$M7,$FX,$crate eq 'zb'?'enclosure; exact cancel gives mid +0 with rad from inputs (PLAN P7)':'exact cancel -> +0 (-0 under Down); near cancel correctly rounded');
  $R->('tie','exact half-ulp ties at p','constructed ties',$hasrm?$X7:'n/a',$FX,$crate eq 'zb'?'enclosure':'ToEven -> even neighbour, ToOdd -> odd neighbour, directed modes per direction, None -> RNDZ@rp(p,64)', '', ($hasrm||$crate eq 'zb'||!$G123{$bare}) ? (r4class=>'tie') : (na=>1)) if $hasrm||$crate eq 'zb';
  if(!$hasrm && $crate eq 'zf'){ if($G123{$bare}){ $R->('tie','exact half-ulp ties at p','constructed ties','n/a','proof','n/a: Spence G1-G3 family (accessor / construct-copy / exact-int-modular): no rounding happens, proven in the committed proof file','',na=>1,g123=>1) }
    else { $R->('tie','exact half-ulp ties at p (family has no RoundingMode parameter: its fixed internal mode)','constructed ties at p in {32,64,128}','fixed mode of the family','fast+heavy',$other?'':'exact tie resolved by the fixed mode the family uses internally (read from the call at '.$src.'), equal to the oracle at that mode; no double rounding','',r4class=>'tie') } }
  # 5: exact range thresholds (EXPONENT_MIN/MAX, not i64)
  my $ovf = $crate eq 'zb' ? 'enclosure, or the all-reals ball (rad=Infinity) on midpoint overflow (CONTRACT 10.4); never a finite ball excluding the true value' : 'ExponentOverflow(sign): ±Inf for modes away from zero; max finite flagged inexact for modes toward zero (ToZero; Down for +; Up for -) per CONTRACT 10.4 description of zf';
  my $unf = $crate eq 'zb' ? 'radius absorbs the lost value and is never zero (CONTRACT 10.4, 3.2)' : 'correctly rounded subnormal / signed zero per mode (EXPONENT_MIN gradual underflow)';
  $R->('ovf-max-finite','result = max finite at p (exponent EXPONENT_MAX = 2^62-1)','operands whose exact result is max finite',$hasrm?$X7:$M7,$FX,$crate eq 'zb'?$ovf:'exact max finite','',r4class=>'huge-exp');
  $R->('ovf-max-plus-half-ulp','exact result = max finite + half ulp (and just below/above)','constructed',$hasrm?$X7:$M7,$FX,$ovf,$crate eq 'zb'?'Ok(all-reals ball)':'ExponentOverflow(Sign)',r4class=>'huge-exp');
  $R->('unf-min-normal','result = smallest normal 2^(EXPONENT_MIN-1)','constructed',$hasrm?$X7:$M7,$FX,$unf,'',r4class=>'tiny-exp');
  $R->('unf-min-subnormal','result = smallest subnormal 2^(EXPONENT_MIN-p)','constructed',$hasrm?$X7:$M7,$FX,$unf,'',r4class=>'subnormal-like');
  $R->('unf-half-min-subnormal-tie','exact result = half the smallest subnormal (tie to 0 or min subnormal)','constructed',$hasrm?$X7:$M7,$FX,$crate eq 'zb'?$unf:'ToEven/ToZero/Down(+) -> +0; Up(+)/FromZero -> min subnormal; ToOdd -> min subnormal; None -> RNDZ (+0)','',r4class=>'tiny-exp');
  if($i64x{$bare}){ $R->('i64-extreme-exponent-arg','exponent argument i64::MIN, i64::MIN+1, i64::MAX-1, i64::MAX, EXPONENT_MAX±1, EXPONENT_MIN±1','constructed',$hasrm?$X7:$M7,$FX,'no wrap; result saturates to ExponentOverflow(sign) / underflow per mode; ilogb/frexp exact','ExponentOverflow(Sign)',keep=>1) }
  # 4: huge trig args
  if($trig{$bare}){ for my $e ('2^16-1','2^16','2^16+1','EXPONENT_MAX'){ my $rv = $crate eq 'zb' ? 'Err(BoundError::Resource { .. }) (zf PrecisionRetryExhausted maps to Resource: zenith-bound/src/mag/convert.rs:209, arith.rs:595, bridge.rs:294)' : $V->('PrecisionRetryExhausted'); my $what = $bare eq 'rem_pi' ? 'x mod pi correctly rounded at p in rm (FIXLIST:31 ZF-17)' : 'correct bits at p';
      $R->("trig-huge-arg e=$e","|x| with exponent $e (trig_reduce_extra cap 2^16, zenith-float-num/src/ops/util.rs:92-95; same cap on x86_64 and i686: 1<<16 with usize saturation)",'x = (1+2^-k)*2^e, random mantissas',$hasrm?$X7:$M7,'fast(ToEven,p=64)+heavy',"$what, or $rv; identical on i686 and x86_64; a wrong value is FAIL (FIXLIST:15 ZF-01, FIXLIST:127 ZF-M22)",$rv,keep=>1) } }
  # 12: T4 mul operand pairs
  if($bare eq 'mul' && $crate eq 'zf'){ for my $w (69,70,71){ my $l=2*$w; $R->("mul-pair-T4 sm=$w lg=$l","unbalanced operands: sm=$w words, lg=$l words (mantissa/mul.rs:53 lg.len()/2 >= sm.len() && sm.len() > 70)","p_a=$w*W, p_b=$l*W, both orders",$hasrm?$X7:$M7,$FX,'correctly rounded at p in rm, equal across the switch','',keep=>1) } }
  # primality (Paul)
  if($bare=~/prime|primality|miller/ && $bare!~/^(ai|bi)_prime$/){ $R->('doc-probable-prime','primality test naming/doc rule (Paul)','docs at '.$src,'n/a','fast','name or docs say "probable prime", or the function is deterministic for all inputs; currently miller_rabin docs (modular.rs:105-108) say neither -> FAIL until FIXLIST doc fix','',keep=>1) }
  # 2: precision axis A to Z
  if($hasp){
    $R->('prec-1-2-3','p in {1,2,3}','normal operands',$hasrm?'all 7 modes':'n/a','fast+heavy',$crate eq 'zb'?'enclosure at prec':$cr,'',r4class=>'prec-boundary');
    $R->('prec-axis','every axis point (i686 96, x86_64 49) incl. k*W-1,k*W,k*W+1 and thresholds T-1,T,T+1','normal operands',$M7,$FX,$crate eq 'zb'?'enclosure at prec':$cr,'',r4class=>'prec-boundary');
    my $rq = $crate eq 'zb' ? 'enclosure or a prompt '.$V->('Resource') : 'correct result, or a prompt refusal '.$V->('MemoryAllocation').' or '.$V->('PrecisionRetryExhausted');
    my $rvv = $crate eq 'zb' ? $V->('Resource') : $V->('MemoryAllocation').' | '.$V->('PrecisionRetryExhausted');
    $R->('prec-MAX-x86_64','p = MAX_PRECISION = 2^60-1 on x86_64','one normal operand','ToEven','heavy (shard pm-max-x86_64; 10 min per row; timeout = FAIL)',"$rq within 10 minutes; a timeout or a wrong value is FAIL and becomes a FIXLIST row. $SKIP",$rvv,keep=>1,vs=>$SKVS);
    $R->('prec-MAX-i686','p = MAX_PRECISION = 2^30-1 on i686 (2^30 bits is allocatable, Paul ruling)','one normal operand','ToEven','heavy only (own shard pm-max-i686; budget 10 min per row; timeout = FAIL)',"$rq within the 10-minute row budget; a timeout or a wrong value is FAIL and every FAIL becomes a FIXLIST row. $SKIP",$rvv,keep=>1,vs=>$SKVS);
    $R->('prec-0-refuse','p = 0','normal operand','ToEven','fast',"refusal: ".$V->('InvalidArgument'),$V->('InvalidArgument'),keep=>1);
    $R->('prec-MAXp1-refuse','p = MAX_PRECISION + 1','normal operand','ToEven','fast',"refusal: ".$V->('InvalidArgument'),$V->('InvalidArgument'),keep=>1);
  } else { $R->('prec-na','signature takes no precision','n/a','n/a','proof','n/a by type: no precision parameter','',na=>1) }
  # 1: seven modes
  for my $m (@MODES){ if($hasrm){ $R->("rm-$m","rounding mode $m vs MPFR $mpfrmode{$m}",'normal + tie + cancellation operands','p=64 fast; axis heavy',$FX,$m eq 'None'?'R6: RNDZ of the exact value at rp(p,64)':$m eq 'ToOdd'?'nearest, exact ties to the odd neighbour':"correctly rounded, $m",'',r4class=>($m eq 'ToOdd'?'rm-ToEven':"rm-$m")) }
    else { $R->("rm-$m","rounding mode $m",'n/a','n/a','proof','n/a by type: no RoundingMode parameter','',na=>1) } }
}
open my $o,'>',$OUT or die; my @H=qw(crate family source category group class class_def inputs modes tier required_result refusal_variant oracle_x86_64 oracle_i686 check verdict_space proof);
sub csvq { my $s=shift; $s="" unless defined $s; $s=~s/\s*\n\s*/ /g; $s=~s/"/""/g; "\"$s\"" }
print $o "# precision-matrix-edges v4 (Harper 2026-10-30). Lines starting with # are metadata; the CSV header is the first line without #.\n# source commit $SRC_COMMIT (generator refuses if zenith-float-num/{src,tests} or zenith-bound/{src,tests} differ)\n";
print $o "# input docs/test-plan/inputs/$_ sha256:$PIN{$_}\n" for sort keys %PIN;
print $o join(',',map {csvq($_)} @H),"\n"; print $o join(',',map {csvq($_)} @$_),"\n" for @out; close $o;
my %fc; $fc{$$_[0]}++ for @fam; my %wp=(); $wp{$$_[1]}=1 for grep {$$_[6]=~/working-precision/} @fam;
print "families zf=$fc{zf} zb=$fc{zb} rows=",scalar(@out),"\n"; print "$_ $cnt{$_}\n" for sort keys %cnt; print "wp families: ",join(' ',sort keys %wp),"\n";
my %oth; $oth{$$_[2]}++ for grep {$$_[0] eq 'zf' && $$_[5] eq 'other'} @fam; print "other keys ",scalar(keys %oth)," bare; r4 fams ",scalar(keys %r4fam),"; missing: ",join(' ',grep {!$r4fam{$_}} sort keys %oth),"\n";

{ my %sh; for my $r (@out){ my ($cr,$fm,$tier,$vs)=@{$r}[0,1,9,15]; my $k = $vs eq 'NA_PROOF' ? 'proof' : $tier=~/pm-max-i686/ ? 'pm-max-i686' : $tier=~/pm-max-x86_64/ ? 'pm-max-x86_64' : $tier=~/fast/ ? 'fast' : 'heavy'; $sh{"$k|$cr|$fm"}++ }
  open my $so,'>',$SHOUT or die; print $so "shard,crate,family,rows\n"; for (sort keys %sh){ my ($a,$b,$c)=split /\|/; print $so join(',',map {csqq($_)} $a,$b,$c,$sh{$_}),"\n" } close $so; print "shard lines: ",scalar(keys %sh),"\n" }
sub csqq { my $s=shift; $s=~s/"/""/g; "\"$s\"" }
open my $n5,'>',$N5OUT or die; print $n5 "crate,family,N,citations\n"; my %fam5;
for my $F (@fam){ my ($crate,$key)=@$F; my $c=$n5{$key}; print $n5 join(',',map {csvq($_)} $crate,$key,$c?5:3,$c?join(' | ',@$c):'none (N=3)'),"\n"; $fam5{$c?5:3}++ } close $n5;
print "N5 families: $fam5{5}; N3 families: $fam5{3}\n";
for my $i (0..$#N5){ my @miss=grep {$_ ne '*' && !$n5hit{$i}{$_}} @{$N5[$i][2]}; print "N5 unmatched names (entry $i): @miss\n" if @miss }
for my $p (@PROOFS){ my ($pp,$key,$src,$bare)=@$p; my $f="$T/$pp"; next if -e $f; open my $pf,'>',$f or die "$f"; print $pf "# NA_PROOF: $key tie row\n\nClaim: `$key` ($src) takes no RoundingMode and performs no rounding (Spence G1-G3: accessor, construct-copy or exact integer/modular), so an exact half-ulp tie at p cannot arise.\n\nEvidence: TO BE WRITTEN by the proving agent (cite the signature and body at $src).\n\nSecond-agent sign-off: UNSIGNED\n\n<!-- pm-coverage FAILs this row until the line above reads `Second-agent sign-off: <agent> <date> <commit>` -->\n"; close $pf }
print "proof placeholders: ",scalar(@PROOFS),"\n";

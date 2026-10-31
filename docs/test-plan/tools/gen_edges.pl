#!/usr/bin/perl
# Harper 2026-10-30: precision-matrix edge list v6 (Spence pmrev5 items 1-14 on top of v5: one hedge rule in tools/pm_hedge.pl, per-class pins, MPFR or mpmath primary oracle).
use strict; use warnings; use FindBin; require "$FindBin::Bin/pm_hedge.pl"; require "$FindBin::Bin/pm_target.pl"; require "$FindBin::Bin/pm_maxp.pl";
my $T=$ENV{PM_TREE}//"$ENV{HOME}/harper-work/zf-precision-matrix";   # worktree; sources pinned to $SRC_COMMIT
my $IN="$T/docs/test-plan/inputs";
my $SRC_COMMIT="7955836a845bd6419dc41d21eda6e739e3b622c5";
my %PIN=("zf-inventory.csv"=>"ba79198eb3486944e3b43a898389b81a7ce8fb69e4f67450e49841e8fdac1d59","zb-inventory.csv"=>"b65c2bad80431500ab6280eb6ca545194407a35e9ae2f9719b5a6e5a38bdcd65","r4-required-results.csv"=>"2746dc4d05cca03d07db6cd8db600e86369a7f33b450c6b679b1d0c1bf0e8c6d","spence-gen_r4.py"=>"aa4d614617d086a49a2d625adbef1b516a7cac7ea37de58e391d40df2b51243f","FIXLIST.csv"=>"8e3a564b5f8cee8782e5d98261aa03f7369dbda35763ffa791592ab8addd07f2","M32-DECISION.md"=>"e6bcf204eaf5f0f57bee9ada064f47a09c66eee658bf3d35c3e343950fb65c5a","td-retriage.csv"=>"fda2e99bf461b128f5145841afcc6d813fdcc79f4d20ce3186efa29e1385b88a","targets.csv"=>"29bc730c658b0d9b6b4d83d5f04f9053161d54c25e6ab2a36ce5afdb7bced09c","tie-resolution.tsv"=>"9d782bf8374727e2fc7c6e2fdcfcb480372f79360ef31b78d33e85bf2519edc4","pm-constants.csv"=>"e77c11f62cc2e656da16b872a01ae8326e315067a0cb45c679d2486191623bdd","fixes_after_1.0.12.txt"=>"03085a0acb3d5d0e51d8e8df13b3c4ad30d7843152d3d32024d6899c9f530c82","zf-audit-SUMMARY.md"=>"a7af866c6b2c6de32f4ed243170a4aa92fff19bc7e3978e86642a3dff1cb75f0");
for my $f (sort keys %PIN){ my $o=`sha256sum "$IN/$f"`; my ($h)=$o=~/^([0-9a-f]{64})/; die "PIN MISMATCH $f: got ".($h//"missing")." want $PIN{$f}\n" unless ($h//"") eq $PIN{$f} }
# Source pin (Harper A2, replaces the A pathspec exclusion): every file under the four source dirs must equal $SRC_COMMIT, except files listed in %SRCPIN, which must equal their pinned sha256.
# zenith-float-num/tests/s9_specials_golds.rs: ZF-T02 (FIXLIST:208) emulated-runner budget scaling (ZF_TIME_FACTOR, integer milli-factor) added in 683f4bb and A2; gen reads no row data from it; any other edit FAILs the pin.
my %SRCPIN=("zenith-float-num/tests/s9_specials_golds.rs"=>"1a30b1fb103272627d2fd8f59b6e6e36ed6da20db0e0e64f261af347867e48f5");
{ my @SD=("zenith-float-num/src","zenith-float-num/tests","zenith-bound/src","zenith-bound/tests");
  my @d=split /\n/, `git -C "$T" diff --name-only $SRC_COMMIT -- @SD`; die "PIN MISMATCH: git diff failed\n" if $?;
  push @d, split /\n/, `git -C "$T" ls-files --others --exclude-standard -- @SD`; die "PIN MISMATCH: git ls-files failed\n" if $?;
  for my $p (@d){ die "PIN MISMATCH: $p differs from $SRC_COMMIT and is not in %SRCPIN\n" unless exists $SRCPIN{$p} }
  for my $p (sort keys %SRCPIN){ my ($h)=`sha256sum "$T/$p"`=~/^([0-9a-f]{64})/; die "PIN MISMATCH $p: got ".($h//"missing")." want $SRCPIN{$p}\n" unless ($h//"") eq $SRCPIN{$p} } }
my %TIE; { open my $th,'<',"$IN/tie-resolution.tsv" or die "tie-resolution.tsv\n"; <$th>; while(<$th>){ chomp; my @c=split /\t/,$_,-1; die "TIE DUP $c[0]\n" if $TIE{$c[0]}; $TIE{$c[0]}=[@c[2..5]] } }
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
# pmrev4 item 3: Paul's corrected FIXLIST rows ZF-WP01..WP10 and ZF-D01, cited by line (generator refuses if absent)
my %WPL; my $D01; { open my $fh,'<',"$T/$FX" or die "$FX"; my $l=0; while(<$fh>){ $l++; $WPL{$2}="$FX:$l ($1)" if /^(ZF-WP\d\d),[^,]*,[^,]*,ExactNum::(\w+):/; $D01="$FX:$l (ZF-D01, LIBRARY.md:399-410 -0i side of the cut)" if /^ZF-D01,/ } }
die "PIN MISMATCH: FIXLIST lacks ZF-D01\n" unless $D01; die "PIN MISMATCH: FIXLIST has ".scalar(keys %WPL)." ZF-WP rows, want 10\n" unless keys %WPL==10;
# pmrev7 (Jeff ruling 2026-10-30 6:00 PM CT): ZF-M27, one gold per row for both targets; cited by line (generator refuses if absent)
my $M27; { open my $fh,'<',"$T/$FX" or die "$FX"; my $l=0; while(<$fh>){ $l++; $M27="FIXLIST:$l ZF-M27" if /^ZF-M27,/ } } die "PIN MISMATCH: FIXLIST lacks ZF-M27\n" unless $M27;
# pmrev8 (FIXLIST 0ef77110, pop-logs/PRINCIPLES.md P1): ZF-M28/M29/M30/D04 cited by line (generator refuses if absent)
my %ZFX; { open my $fh,'<',"$T/$FX" or die "$FX"; my $l=0; while(<$fh>){ $l++; $ZFX{$1}="FIXLIST:$l $1" if /^((?:ZF-(?:M28|M29|M30|M31|M32|D04|D07|D08|G01)|ZB-D0[56])),/ } } for (qw(ZF-M28 ZF-M29 ZF-M30 ZF-M31 ZF-M32 ZF-D04 ZF-D07 ZF-D08 ZF-G01 ZB-D05 ZB-D06)){ die "PIN MISMATCH: FIXLIST lacks $_" unless $ZFX{$_} } my ($M28,$M29,$M30,$D04)=@ZFX{qw(ZF-M28 ZF-M29 ZF-M30 ZF-D04)};
my @TWIN; my %M30DONE; my $M30RE=qr/^(?:ExactNum|ExactComplex|ExactNumArray|Ball|ComplexBall|\$arr\[ieee_soft\/array\.rs\])::(?:sin|cos|tan|atan|asin|acos|ln|log|exp|sinh|cosh|tanh|pow|atan2|log2|log10|sin_cos)$/;
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
 (map { my ($n,$s,$d)=@$_; [q(zf),q(.*),[$n],"$FX:18 (ZF-04 S1 wrong-answer, 29 = 25 catalog.rs + 4 heun_fox.rs, Paul corrected 2026-10-30); zenith-float-num/src/$s at SRC_COMMIT 7955836 (doc: Accuracy: unspecified)"] } ([qw(scorer_gi catalog.rs:575-576)],[qw(scorer_hi catalog.rs:603-604)],[qw(kelvin_ber_nu catalog.rs:629-630)],[qw(kelvin_bei_nu catalog.rs:658-659)],[qw(kelvin_ker_nu catalog.rs:687-688)],[qw(kelvin_kei_nu catalog.rs:722-723)],[qw(kelvin_ber catalog.rs:757-758)],[qw(kelvin_bei catalog.rs:772-773)],[qw(kelvin_ker catalog.rs:787-788)],[qw(kelvin_kei catalog.rs:802-803)],[qw(struve_h catalog.rs:821-822)],[qw(anger_j catalog.rs:916-917)],[qw(weber_e catalog.rs:962-963)],[qw(clausen_cl2 catalog.rs:1007-1008)],[qw(clausen_cl3 catalog.rs:1038-1039)],[qw(hurwitz_zeta catalog.rs:1066-1067)],[qw(riemann_zeta catalog.rs:1115-1116)],[qw(polygamma catalog.rs:1148-1149)],[qw(barnes_g catalog.rs:1181-1182)],[qw(hypergeom_0f1 catalog.rs:1328-1329)],[qw(hypergeom_1f1 catalog.rs:1347-1348)],[qw(hypergeom_pfq catalog.rs:1389-1390)],[qw(lambert_w0 catalog.rs:1469-1470)],[qw(lambert_wm1 catalog.rs:1484-1485)],[qw(polylog catalog.rs:1581-1582)],[qw(meijer_g heun_fox.rs:576-577)],[qw(fox_h heun_fox.rs:640-641)],[qw(heun_g heun_fox.rs:771-772)],[qw(heun_c heun_fox.rs:967-968)])), (map { [q(zf),q(.*),[$_->[0]],"$FX:18 (ZF-04 range catalog.rs:575-1581); zenith-float-num/src/$_->[1] at SRC_COMMIT 7955836 (doc: Accuracy: correctly-rounded, NOT unspecified; outside the 29 unspecified, kept at N=5 because it is in the cited range and goes through round_exact)"] } ([qw(jacobi_arcsn catalog.rs:1254-1255)],[qw(jacobi_arccn catalog.rs:1273-1274)],[qw(jacobi_arcdn catalog.rs:1300-1301)])),
 (map { [q(zf),q(.*),[$_],"$WPL{$_} (S1 wrong-answer: only 32 bits right at every p)"] } sort keys %WPL),
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
my $UNP='[PM1-UNPINNED] required result: FAIL until the gold pins one outcome per row at PM-1; pin rule: ';
my $PIB='the documented up-front pi-bit budget constant (ZF-M22, FIXLIST:127; set by Rusty, the same value on x86_64 and i686)';
my $OVFZ='ExponentOverflow(sign): ±Inf for modes away from zero; max finite flagged inexact for modes toward zero (ToZero; Down for +; Up for -) per CONTRACT 10.4 description of zf';
my %HXEX=map{$_=>1} qw(abs neg copysign floor ceil round trunc fract clamp max min signum sign);
my %HXCR=map{$_=>1} qw(atan atan2 tanh erf erfc sqrt cbrt ln log log2 log10 log1p asinh acosh normal_cdf ln_gamma lambert_w0);
my %HXINV=map{$_=>1} qw(asin acos atanh);
my %ZPOLE=(kelvin_ker=>1); my %ZUNP=map{$_=>1} qw(kelvin_ker_nu kelvin_kei_nu fox_h meijer_g);
my %VALCL=map{$_=>1} qw(normal negative zero signed-zero cancellation tie tiny-exp subnormal-like unf-min-normal unf-min-subnormal unf-half-min-subnormal-tie prec-1-2-3 prec-axis prec-MAX-i686);
sub Vk { my ($a,$v)=@_; return "Err(BoundError::$v)" if $$a{crate} eq 'zb'; my $r=$$a{_ret}//'value'; return 'None' if $r eq 'option'; return "Err(Error::$v)" if $r eq 'result'; return "false (documented refusal, $v-class input)" if $r eq 'bool'; "NaN whose error kind is Error::$v" }
my %i64x=map {$_=>1} qw(scalb ldexp ilogb frexp);
my $W=0;
my @MODES=qw(None Up Down ToZero FromZero ToEven ToOdd);
my %mpfrmode=(None=>'RNDZ at rp(p,64) (R6)',Up=>'RNDU',Down=>'RNDD',ToZero=>'RNDZ',FromZero=>'RNDA',ToEven=>'RNDN',ToOdd=>'RNDN, exact tie flipped to the odd neighbour by the oracle (tie detected as exact at p+1 bits)');
my @out; my %cnt; our @PROOFS;
# pmrev4 item 4: one outcome per executed row. Every either/or required result is rewritten to a single pinned
# outcome (one value or one named refusal; the other branch is FAIL) or to a [PM1-UNPINNED] pin rule that FAILs
# in pm-coverage rule d until the gold pins it. The generator refuses to write a table with either/or text left.
our @LEFT;
sub pin1 { my $a=shift; return if ($$a{verdict_space}//'') eq 'NA_PROOF';
  my ($k,$c,$cr)=@$a{qw(family class crate)}; my $o=$$a{oracle_x86_64}//''; my $bare=$k=~s/.*:://r;
  local $_=$$a{required_result}//''; my $v=$$a{refusal_variant}//'';
  my $U='[PM1-UNPINNED] required result: FAIL until the gold pins one outcome per row at PM-1; pin rule: ';
  my $conv=' [convention: mpmath quad/findroot/odefun; LAPACK getrf/geqrf/gesvd normalisations; Spence pmrev4 item 4 rulings]';
  if(/UNCLEAR/){
    if($k=~/^gauss_kronrod::|^roots::(brent|illinois)$/){ $_='R4 result correctly rounded at p in rm (Ziv-certified against the pinned oracle); a refusal is FAIL'.($k=~/^roots/?' (Error::InvalidArgument is pinned only on the no-sign-change row)':'').$conv; $v='value' }
    elsif($k eq 'ode::rk4'){ $_='R4 the exact fixed-step RK4 recurrence for the given h and step count, every step correctly rounded at p in rm; a refusal is FAIL (Error::InvalidArgument is pinned only on the bad-step-interval row)'.$conv; $v='value' }
    elsif($k eq 'ode::rk45_adaptive'){ $_=$U.'the value within tol of mpmath odefun at 2p bits in env oracle-2026-10-30 where the gold records that the oracle meets tol; NaN whose error kind is Error::PrecisionRetryExhausted where the gold records that it does not'.$conv; $v='pinned per row in the gold' }
    elsif($k=~/_decomp$/){ $_='R4 Some(factors): each factor entry correctly rounded at p of the exact factor under the pinned normalisation (LU: partial pivoting, first max abs(pivot), lowest index on exact ties; QR: diag(R)>=0; SVD: Sigma descending, first nonzero of each u_i positive; eig: ascending, unit-norm, sign-normalised); the generated inputs exclude the degenerate cases, so None is FAIL (None is pinned only on the degenerate row)'.$conv; $v='Some(factors)' }
    else { push @LEFT,"UNCLEAR without a ruling: $k $c"; } }
  s/ or the defined constructor value/ (the stored field is the defined constructor value)/;
  if(/^(R4 unnormalized input either|from_raw_parts normalises exactly)/){ $_='from_raw_parts with these parts (valid at p: exponent in range, limb count matching p) normalises exactly; the family result is then '.($cr eq 'zb'?'an enclosure of the exact value':'correctly rounded at p in rm').'; a refusal is FAIL; a value computed from the unnormalised mantissa is FAIL (invalid parts: raw-parts-invalid row of ExactNum::from_raw_parts, pinned Error::InvalidArgument) [convention: from_raw_parts_checks.rs; Spence pmrev4 item 4]'; $v='' }
  if(/^correct limit \(IEEE 754-2019 7\.2\/9\.2, MPFR semantics\) or refusal$/ || /^R4 IEEE-like limit value/){ my ($cv)=/(\[convention:[^\]]*\])$/; $_='the limit the pinned x86_64 oracle returns for these +-Inf operands at p in rm (IEEE 754-2019 7.2/9.2); where that oracle returns NaN the pin is the refusal in refusal_variant; any other result is FAIL'.($cv?" $cv":''); $v=~s/ where the limit is undefined$/ (pinned where the oracle returns NaN)/ }
  s/^enclosure, or the all-reals ball \(rad=Infinity\) (?:on midpoint overflow|when the midpoint overflows) \(CONTRACT 10\.4\)/enclosure containing the exact value; where the exact midpoint exceeds the max finite value the pinned result is the all-reals ball (rad=Infinity, CONTRACT 10.4)/;
  s/; never a finite ball excluding the true value/; a finite ball excluding the true value is FAIL/;
  if($c eq 'prec-MAX-x86_64'){ my ($sk)=/(MAX-timeout skip .*)$/; my $mv=(split / \| /,$v)[0]; $_="a prompt ".($cr eq 'zb'?'':'refusal ')."$mv within 10 minutes (p = 2^60-1 bits = 2^57 bytes cannot be allocated, so any returned value is FAIL); a timeout is FAIL and becomes a FIXLIST row. $sk"; $v=$mv }
  if($c eq 'prec-MAX-i686'){ my ($sk)=/(MAX-timeout skip .*)$/; $_=($cr eq 'zb'?'enclosure at prec':'correct result at p')." = 2^30-1 (allocatable on i686, Paul ruling) within the 10-minute row budget; a refusal is FAIL, a timeout is FAIL, a wrong value is FAIL, and every FAIL becomes a FIXLIST row. $sk"; $v='' }
  if(/^R3: correctly rounded at the requested p in the requested rm, or refuse with/){ $_='R3: correctly rounded at the requested p in the requested rm (valid operands; invalid operands are the out-of-domain row); a refusal is FAIL; a result accurate only to 32 bits is FAIL ('.($WPL{$bare}//'FIXLIST ZF-WP row MISSING').')'; push @LEFT,"no ZF-WP row for $k" unless $WPL{$bare} }
  s/length 0 or 1 window/length-0 window, length-1 window/;
  s/pollard_rho prime or cap hit -> None/pollard_rho on a prime -> None; pollard_rho cap hit -> None/;
  s/exact 0 \(pdf\/pmf\) or 0\/1 \(cdf\)/exact 0 for pdf\/pmf; cdf: exact 0 below the support, exact 1 above it/;
  if(/^R4 exponent extremes: no wrap/){ my ($cv)=/(\[convention:[^\]]*\])$/; $_=$U."no wrap (checked arithmetic); the value at the cap this family's doc states, recorded per row with the refusal named in refusal_variant where the doc states a refusal $cv" }
  if(/^R4 zero polynomial/){ my ($cv)=/(\[convention:[^\]]*\])$/; $_=$U."degree and derivative of the zero polynomial exactly as the poly.rs doc states them (one representation, recorded in the gold); integral constant 0 $cv" }
  s{oracle MPFR or mpmath Ziv pair p\+256/p\+512}{'oracle '.($o=~/^primary: MPFR/ ? 'MPFR (as pinned in oracle_x86_64)' : $o=~/^primary: mpmath 1\.3\.0 exact integer/ ? 'mpmath exact integers (as pinned in oracle_x86_64)' : $o=~/digest:/ ? 'gold digest (as pinned in oracle_x86_64)' : 'mpmath Ziv pair p+256/p+512 (as pinned in oracle_x86_64)')}e;
  s/values are exact or correctly rounded at p in rm/zero, infinite and listed exact values are exact; every other value is correctly rounded at p in rm/;
  s/-0 gives -Inf \(odd n\) or \+Inf \(even n\)/-0 gives -Inf for odd n and +Inf for even n/;
  s/\(FIXLIST doc row for LIBRARY\.md:399 is Paul's\)/(doc row: $D01)/;
  if($c eq 'doc-probable-prime'){ $_='FAIL until the FIXLIST doc fix lands: the name and docs must say "probable prime", because a Miller-Rabin test is not proven deterministic for all inputs here; today modular.rs:105-108 does not say it' }
  if($v eq 'DivisionByZero or None per doc'){ $_=$U.$_.'; refusal per family: the one the fn doc states for modulus 0 (mod_pow, mod_inv), none for the other cases'; $v='pinned per row from the fn doc' }
  $v=~s/^value or NaN\(DivisionByZero\)$/value at a non-pole zero; NaN whose error kind is Error::DivisionByZero at a pole/;
  $v=~s/^entries or None$/Some(factors)/;
  $$a{required_result}=$_; $$a{refusal_variant}=$v }
sub pin2 { my $a=shift; my ($k,$c,$cr)=@$a{qw(family class crate)}; my $bare=$k=~s/.*:://r;
  local $_=$$a{required_result}//''; my $v=$$a{refusal_variant}//''; my ($cv)=/(\[convention:[^\]]*\])$/; my $CV=$cv?" $cv":'';
  my $ex=($$a{verdict_space}//'') ne 'NA_PROOF'; my $done=0;
  if($ex && $cr eq 'zf' && $k=~/^ExactNumArray::(lu|qr|svd|eigen)_decomp$/){ my $d=$1; my $dc=' [convention: array.rs:1315-1497 (Option return); LAPACK getrf, geqrf, gesvd normalisations; Spence pmrev5 items 6-7]';
    my $norm='(LU: partial pivoting on the first max abs(pivot); exact inputs, so a tie always goes to the lowest index; QR: diag(R)>=0; SVD: Sigma descending, first nonzero of each u_i positive; eigen: ascending, unit-norm, sign-normalised)';
    if($c=~/^(zero|signed-zero)$/){ my %w=(lu=>'the zero matrix has a zero pivot after partial pivoting (exactly singular)',qr=>'qr of the zero matrix has R_ii = 0',svd=>'the zero matrix has only a repeated zero singular value',eigen=>'the zero matrix has only a repeated zero eigenvalue'); $_="None: $w{$d}; Some(factors) is FAIL$dc"; $v='None' }
    elsif($c=~/^(inf|nan)/){ $_="None: an Inf entry and a NaN entry have no exact factorisation; Some(factors) is FAIL$dc"; $v='None' }
    elsif($c eq 'out-of-domain'){ $_="None: a non-square matrix and an empty matrix".($d eq 'lu'?' and an exactly singular matrix (a zero pivot after partial pivoting)':'')." have no factorisation under the pinned normalisation; Some(factors) is FAIL$dc"; $v='None' }
    elsif($c eq 'ovf-max-plus-half-ulp'){ $_=$UNP."the single outcome the gold records from the fn doc for the exact factorisation whose largest entry is max finite + half ulp$dc"; $v='pinned per row in the gold' }
    elsif($c=~/^(degenerate|prec-0-refuse|prec-MAXp1-refuse|prec-MAX-x86_64)$/){ }
    else { $_="Some(factors): each factor entry correctly rounded at p in rm of the exact factor under the pinned normalisation $norm; inputs are generated nonsingular with distinct singular values and eigenvalues, so None is FAIL (None is pinned on the zero, signed-zero, inf, NaN, out-of-domain and degenerate rows)$dc"; $v='Some(factors)' }
    $done=1 }
  elsif($ex && $cr eq 'zf' && $k eq 'ode::rk45_adaptive'){ my $h='rk45 ruling (Spence pmrev4, pmrev5 item 4): within tol of mpmath odefun at 2p bits in oracle-2026-10-30, pinned per row; this row: ';
    my $api=' (rk45_adaptive returns Option at 7955836, ode.rs:439, so this row FAILs until it reports the kind: FIXLIST API row needed)';
    if($VALCL{$c} || $c=~/^rm-/){ $_=$h.'the value within tol of mpmath odefun at 2p bits; a refusal is FAIL'; $v='value' }
    elsif($c=~/^prec-(0|MAXp1)-refuse$/){ $_=$h.'refusal Error::InvalidArgument (p outside [1, MAX_PRECISION])'.$api; $v='Error::InvalidArgument' }
    elsif($c eq 'prec-MAX-x86_64'){ $_=$h.$_ }
    else { $_=$h."refusal Error::PrecisionRetryExhausted ($c inputs: odefun at 2p bits cannot meet tol)".$api; $v='Error::PrecisionRetryExhausted' }
    $_.=' [convention: mpmath odefun; Spence pmrev5 item 4]'; $done=1 }
  elsif($ex && $cr eq 'zf' && $k eq 'ExactNum::from_raw_parts'){ my $h='from_raw_parts ruling (Spence pmrev4, pmrev5 item 5): exact normalized, pinned per row; this row: ';
    if($c=~/^(raw-parts-invalid|out-of-domain|inf|nan|ovf-max-plus-half-ulp)/){ $_=$h.'refusal Err(Error::InvalidArgument): these parts do not describe a finite value at p'.($c eq 'ovf-max-plus-half-ulp'?' (exponent above EXPONENT_MAX)':'').'; any returned value is FAIL'; $v='Err(Error::InvalidArgument)' }
    else { $_=$h.'the exact normalised value of these valid parts (from_raw_parts normalises exactly, no rounding); a refusal is FAIL (invalid parts pin Err(Error::InvalidArgument) on the raw-parts-invalid row)'; $v='' }
    $_.=' [convention: ext.rs:1970 Result return; from_raw_parts_checks.rs]'; $done=1 }
  if($ex && !$done && $cr eq 'zf' && $c eq 'huge-exp' && ($trig{$bare} || /\belse\b/i)){ my $cx=$k=~/^ExactComplex::/;
    if($trig{$bare}){ my $m=Vk($a,'MemoryAllocation'); $_="refusal: $m raised before allocating: reducing an argument with exponent in [EXPONENT_MAX-64, EXPONENT_MAX] needs more pi bits than $PIB; identical on i686 and x86_64; any returned value is FAIL (FIXLIST:127 ZF-M22)"; $v=$m }
    elsif($HXEX{$bare}){ $_='exact result (no rounding; the huge exponent is representable); a refusal is FAIL'.$CV; $v='' }
    elsif($HXCR{$bare} || ($cx && $HXINV{$bare})){ $_='correctly rounded at p in rm (for an in-domain operand of this exponent the result is finite and in range); a refusal is FAIL'.$CV; $v='' }
    elsif($HXINV{$bare}){ my $m=Vk($a,'InvalidArgument'); $_="refusal: $m (|x| > 1 is outside the domain); any returned value is FAIL$CV"; $v=$m }
    else { $_=$UNP.'the one outcome the gold records for this operand: the correctly rounded value at p in rm where the exact result is in range, '.Vk($a,'ExponentOverflow').' with the sign of the exact result where it exceeds EXPONENT_MAX'.$CV; $v='pinned per row in the gold' }
    $done=1 }
  if($ex && !$done && $cr eq 'zf' && $c=~/^ovf-/ && /\belse\b/i){ if($c eq 'ovf-max-finite'){ $_='R4 exact max finite value at p (representable, no rounding); a refusal is FAIL'.$CV; $v='' } else { $_="R4 $OVFZ$CV"; $v='ExponentOverflow(Sign)' } $done=1 }
  if($ex && !$done && $cr eq 'zf' && $c=~/^(tiny-exp|unf-min-normal|unf-min-subnormal|unf-half-min-subnormal-tie|subnormal-like)$/ && (/\belse\b/i || m{ / } || $v=~m{/})){
    $_ = $c eq 'unf-min-normal' ? 'exact smallest normal 2^(EXPONENT_MIN-1) (representable, no rounding); a refusal is FAIL' : $c eq 'unf-half-min-subnormal-tie' ? 'the correctly rounded value in rm of half the smallest subnormal (ToEven, ToZero, Down(+) give +0; Up(+), FromZero, ToOdd give the min subnormal; None gives +0)' : $c eq 'unf-min-subnormal' ? 'exact smallest subnormal 2^(EXPONENT_MIN-p) (gradual underflow at EXPONENT_MIN, representable); a refusal is FAIL' : 'the correctly rounded value at p in rm under gradual underflow at EXPONENT_MIN (a subnormal result and a signed-zero result are both that rounded value); a refusal is FAIL';
    $_.=$CV unless index($_,'[convention:')>=0; $v=''; $done=1 }
  if($ex && !$done && $cr eq 'zf' && $c eq 'out-of-domain'){ my @t=grep {length} map {s/^\s+|\s+$//gr} split /;/,$v;
    if(@t>1 || /\belse\b/i || m{cap/budget}){ my $r1=(grep {/InvalidArgument/} @t)[0] // Vk($a,'InvalidArgument'); $_="refusal $r1 for inputs outside the documented domain of this family (shape, index, radix, p and real-domain violations; exact poles are not inputs of this row); any returned value is FAIL$CV"; $v=$r1; $done=1 } }
  if($ex && !$done && $v=~/non-pole zero/){
    if($ZPOLE{$bare}){ my $m=Vk($a,'DivisionByZero'); $_="the pole at 0: +Inf with error kind Error::DivisionByZero (log singularity, DLMF 10.61); any finite value is FAIL$CV"; $v=$m }
    elsif($ZUNP{$bare}){ $_=$UNP."the one outcome the gold records at 0 for the given parameters: the limit value where the family is finite at 0, NaN with Error::DivisionByZero where it has a pole at 0$CV"; $v='pinned per row in the gold' }
    else { $_="the exact limit value at 0 (DLMF), correctly rounded at p in rm; this family has no pole at 0, so a refusal is FAIL$CV"; $v='' } $done=1 }
  for my $t ($_,$v){
    $t=~s{\bn/a\b}{not applicable}g; $t=~s{IEEE/MPFR}{IEEE and MPFR}g; $t=~s{i686/x86_64}{i686 and x86_64}g; $t=~s{Inf/Inf}{Inf divided by Inf}g;
    $t=~s{compare/contains with NaN -> unordered/false}{compare with NaN -> unordered; contains with NaN -> false}g;
    $t=~s{p=0 / MAX\+1 refused}{p=0 and p=MAX+1 are refused on the prec-0-refuse and prec-MAXp1-refuse rows}g;
    $t=~s{MPFR doc / mpmath Ziv}{MPFR doc, mpmath Ziv}g;
    $t=~s{gmpy2 / python-flint|python-flint/gmpy2}{mpmath exact integers (gmpy2 backend), python-flint second check}g;
    $t=~s{ / }{, }g if $t=~/^(empty input|R4 empty input)/;
    $t=~s{mode-independent where no rm; else values rounded per mode}{$$a{_hasrm}?'values rounded per mode':'mode-independent values'}ge;
    $t=~s{\(-0 only in Down\), else correctly rounded}{(-0 only in Down); a non-cancelling pair gives the correctly rounded value}g;
    $t=~s{; else signed Inf}{; every other Inf operand pair gives the signed Inf}g;
    $t=~s{ / }{, }g;
    $t=~s{(\S+)}{ my $w=$1; (index($w,'/')<0 || pm_slash_ok($w)) ? $w : ($w=~s{(?<=[A-Za-z0-9_)+\-])/(?=[A-Za-z(+\-])}{, }gr) }ge; }
  $v='' if $v=~m{^(value , signed 0|value / signed 0|subnormal, \+-0|subnormal/\+-0)$};
  $v='signed zero (exact)' if $v=~m{^-0 [/,] \+0$};
  if($v=~m{^ExactNum::nan\(.*\)(?: /|,) Err\(}){ $v = $c=~/refuse/ ? Vk($a,'InvalidArgument') : '' }
  $v='None' if $v=~/^None(?:; None)+$/;
  $$a{required_result}=$_; $$a{refusal_variant}=$v }
# pmrev7 (Jeff ruling 2026-10-30 6:00 PM CT, ZF-M27): one gold per row for i686 and x86_64; MAX split resolved per design s4; M27 pin note.
my @F17=qw(crate family source category group class class_def inputs modes tier required_result refusal_variant oracle_x86_64 oracle_i686 check verdict_space proof);
sub pin3 { my $a=shift; my $c=$$a{class}; my $ex=($$a{verdict_space}//'') ne 'NA_PROOF'; my $cr=$$a{crate}; my $zb=$cr eq 'zb';
  my $SK='MAX-timeout skip (Paul, Spence pmrev3 6; one row, both targets): the first target run of this row in a family runs first; if it times out, the run on the other target is not started and is recorded SKIPPED_AFTER_TIMEOUT, which counts as FAIL in the verdict space and in pm-coverage and is never omitted; the family\'s RC stays blocked until its FIXLIST fix lands, and every skipped run then reruns in full';
  my $TIER='heavy only (own shard pm-max, run on both targets; budget 10 min per run; timeout = FAIL)';
  if($c eq 'prec-MAX-i686'){ $$a{class}='prec-2^30-1'; $$a{class_def}='p = 2^30-1 on both targets (MAX_PRECISION on every target after ZF-M28, P1; 2^30 bits is allocatable, Paul ruling)'; $$a{tier}=$TIER;
    local $_=$$a{required_result}; s/ = 2\^30-1 \(allocatable on i686, Paul ruling\)/ = 2^30-1, one gold for both targets ($M27; $M28: MAX_PRECISION = 2^30-1 on every target, so this is the shared limit and must give the correct result)/; s/MAX-timeout skip .*$/$SK/ or $_.=" $SK"; $_.=" ($M28: MAX_PRECISION = 2^30-1 on every target, so p = 2^30-1 is the shared limit and must give this result.)" unless /ZF-M28/; $$a{required_result}=$_ }
  if($c eq 'prec-MAX-x86_64'){ $$a{class}='prec-2^60-1'; $$a{class_def}='p = 2^60-1 on both targets (above MAX_PRECISION = 2^30-1 after ZF-M28, P1 one limit)'; $$a{tier}=$TIER;
    my ($ia,$ma)=$zb?('BoundError::InvalidArgument','BoundError::Resource'):('Error::InvalidArgument','Error::MemoryAllocation');
    $$a{required_result}="refusal: $ia (p = 2^60-1 is above MAX_PRECISION = 2^30-1, one limit on every target after $M28, P1); any returned value is FAIL and $ma is FAIL. Red until ZF-M28 lands: at 7955836 MAX_PRECISION is 2^60-1 on some targets (defs.rs:77".($zb?'; zb mirrors it at zenith-bound/src/precision.rs:16':'')."). $SK"; $$a{refusal_variant}=$ia;
    push @TWIN, {%$a, class=>'prec-2^30', class_def=>'p = 2^30 on both targets (MAX_PRECISION + 1 after ZF-M28, P1 one limit)', required_result=>"refusal: $ia (p = 2^30 = MAX_PRECISION + 1 on every target after $M28, P1); any returned value is FAIL and $ma is FAIL. Red until ZF-M28 lands: at 7955836 2^30 is accepted on some targets (defs.rs:77".($zb?'; zb precision.rs:16':'')."). $SK".($$a{family}=~/::(?:brent|illinois|rk4|rk45_adaptive)$/ ? " Red until $M29 lands: this family returns Option at 7955836; ZF-M29 makes it return Result with the pinned error kind." : '')} }
  if($$a{family}=~/::(?:brent|illinois|rk4|rk45_adaptive)$/){ local $_=$$a{required_result};
    s/ \(rk45_adaptive returns Option at 7955836, ode\.rs:439, so this row FAILs until it reports the kind: FIXLIST API row needed\)/ (red until $M29 lands: rk45_adaptive returns Option at 7955836, ode.rs:439; ZF-M29 makes it return Result with this kind)/g;
    s/; rk4 returns Option at 7955836 \(ode\.rs:73\), so this row FAILs until rk4 reports the kind: FIXLIST API row needed/; red until $M29 lands: rk4 returns Option at 7955836 (ode.rs:73); ZF-M29 makes it return Result with this kind/g;
    if(s/NaN whose error kind is (Error::\w+), held in the returned Option \(roots\.rs:238, 407 return Option<ExactNum>\)/$1 (red until $M29 lands: roots.rs:238, 407 return Option<ExactNum> at 7955836; ZF-M29 makes them return Result with this kind)/g){ $$a{refusal_variant}=~s/^NaN whose error kind is //; }
    $_.=" Red until $M29 lands: this family returns Option at 7955836; ZF-M29 makes it return Result with the pinned error kind." if $ex && $$a{refusal_variant}=~/Error::/ && !/ZF-M29/;
    $$a{required_result}=$_ }
  push @LEFT,"API row $$a{family} $$a{class}" if "$$a{required_result} $$a{refusal_variant}"=~/API row/;
  $$a{class_def}=~s/every axis point \(i686 96, x86_64 49\)/every axis point of the union of the i686 (96) and x86_64 (49) axis sets, each run on both targets with one gold/ if $c eq 'prec-axis';
  $$a{required_result}=~s/the pinned x86_64 oracle/the pinned oracle/;
  if($ex){ $$a{oracle_i686}="same gold as oracle_x86_64: one gold per row for both targets (Jeff ruling 2026-10-30 6:00 PM CT, $M27); the i686 run must be bit-identical to it".($zb?' (byte-identical BallRow digest)':'') }
  if($ex && pm_m27_candidate(@$a{@F17})){ $$a{required_result}.=" Pin note ($M27, Jeff ruling 2026-10-30 6:00 PM CT): cells in the 6 directed modes (ToEven, Up, Down, ToZero, FromZero and ToOdd: every mode except None, which is target-identical per Jeff 2026-10-29) at a p that is not a multiple of 64 are red (FAIL) until ZF-M27 lands; then every mode rounds at rp(p,64) and the result is bit-identical on i686 and x86_64." } }
# pmrev10 (Jeff 2026-10-30 6:31 PM CT, PRINCIPLES.md P6 sha256 d39b9a60, FIXLIST ZF-M31/M32): 64-bit targets only.
# Rows that exist only because of i686/32-bit are dropped (counted in %DROP); every other row is rewritten to name no 32-bit target.
our %DROP; my $P6='PRINCIPLES.md P6 (Jeff 2026-10-30 6:31 PM CT)';
sub pin4 { my $a=shift; my $c=$$a{class}; my $ex=($$a{verdict_space}//'') ne 'NA_PROOF'; my $zb=$$a{crate} eq 'zb';
  if($c eq 'prec-2^30-1'){ $DROP{'prec-2^30-1: i686 MAX_PRECISION (was prec-MAX-i686); ZF-M28 2^30-1 superseded by ZF-M32'}++; return 0 }
  if($c eq 'prec-2^30'){ $DROP{'prec-2^30: ZF-M28 MAX+1 twin; superseded by ZF-M32 (prec-MAXp1-refuse covers MAX+1)'}++; return 0 }
  if($c eq 'm30-p2^24'){ $DROP{'m30-p2^24: i686 calc_mul_cost as u32 threshold (FIXLIST:148); no 64-bit edge'}++; return 0 }
  if($c eq 'prec-2^60-1'){ my $ia=$zb?'Err(BoundError::InvalidArgument)':'Error::InvalidArgument'; $$a{class}='prec-over-MAX';
    $$a{class_def}="p = 2*PM_MAX_P + 1, over-MAX: derived from PM_MAX_P ($main::PMNOTE)";
    $$a{required_result}="refusal: $ia raised up front by the precision check (zf p_assertion num.rs:55-57 against MAX_PRECISION; zb precision.rs:17), before any allocation; any returned value is FAIL; a timeout is FAIL; a crash is FAIL; an OOM kill is FAIL".($zb?" (zb Precision is usize, precision.rs:9; $ZFX{'ZB-D05'} fixes the CONTRACT u32 text)":'').". Red until $ZFX{'ZF-M32'} lands: the pinned source predates M32, so the p this row derives may be accepted and reach allocation. MAX-timeout: a run past the 10-minute budget is recorded SKIPPED_AFTER_TIMEOUT, which counts as FAIL in the verdict space and in pm-coverage and is never omitted";
    $$a{refusal_variant}=$ia }
  if($c eq 'm30-p2^27'){ $$a{class_def}="p = 2^27 (ZF-M30 64-bit cost-overflow edge: series.rs:156,238,244 and tan.rs:128 iter_cost*niter panics at p>=2^27 on the host, $ZFX{'ZF-M30'}; an input point; P1 one gold)" }
  if($c eq 'prec-axis'){ $$a{class_def}="every axis point of the 64-bit axis set (49 points) incl. k*W-1,k*W,k*W+1 with W=64 and thresholds T-1,T,T+1; the W=32 word points are dropped ($P6)" }
  for my $k (keys %$a){ next if index($k,"_")==0; next unless defined $$a{$k}; local $_=$$a{$k};
    s/ Pin note \(FIXLIST:\d+ ZF-M27, Jeff ruling[^)]*\): cells in the 6 directed modes .*?bit-identical on i686 and x86_64\./ ZF-M27 ($M27) is met on every 64-bit target: the word is 64 bits, so every mode already rounds at rp(p,64) ($P6)./;
    s/MAX-timeout skip \(Paul, Spence pmrev3 6; one row, both targets\):.*?reruns in full/MAX-timeout: a run past the 10-minute budget is recorded SKIPPED_AFTER_TIMEOUT, which counts as FAIL in the verdict space and in pm-coverage and is never omitted; the family's RC stays blocked until its FIXLIST fix lands/;
    s/MAX-timeout skip \(Paul, Spence pmrev3 6\):.*?(?:reruns in full|$)/MAX-timeout: a run past the 10-minute budget is recorded SKIPPED_AFTER_TIMEOUT, which counts as FAIL and is never omitted/;
    s/ \(ZF-M19 add panics on i686\)/ (ZF-M19)/g;
    s/ Red until ZF-M27 lands: the p=32 cell is FAIL until rounding happens at rp\(p,64\)\./ The p=32 cell is not red: on every 64-bit target the word is 64 bits, so rounding already happens at rp(p,64) ($M27, $P6)./g;
    s/the same value on x86_64 and i686/one value on every 64-bit target/g;
    s/\b(?:i686 and x86_64|x86_64 and i686)\b/every 64-bit target/g;
    s/\b(?:on |for )?both targets\b/on every 64-bit target/g;
    $$a{$k}=$_ }
  if($ex){ $$a{oracle_i686}="one committed gold (the oracle column) for every active target in docs/test-plan/inputs/targets.csv ($P6, P1); each active target's run must be bit-identical to it".($zb?' (byte-identical BallRow digest)':'')."; a pending target is reported NOT-RUN, never PASS" }
  else { $$a{oracle_i686}='n/a (proof row)' }
  1 }

# pmrev12 (Paul M32 ruling 2026-10-30 7:18 PM CT): MAX, MAX+1 and over-MAX rows reference the named pin PM_MAX_P (inputs/pm-constants.csv, owner Rusty), never a literal.
# Basic ops (add, sub, mul, div, sqrt, from_*, to_*, parse, format, fmt; categories conversion and io/format) pin the correct result at PM_MAX_P.
# Every other family's PM_MAX_P row is TRACE-PENDING (RC-FAIL) until the cross-file allocation trace (design s28.4) pins its outcome.
our $PMX=pm_maxp($T,$SRC_COMMIT,"$IN/pm-constants.csv");
our $PMNOTE="PM_MAX_P is zf MAX_PRECISION, read by tools/pm_maxp.pl from $$PMX{file}:$$PMX{line} at the source commit (value in the CSV header; expected M32 value from inputs/pm-constants.csv, owner $$PMX{owner})";
sub m32basic { my ($fam,$cat)=@_; my ($bare)=$fam=~/([^:]+)$/; return ($bare=~/^(?:add|sub|mul|div|sqrt|parse|format|fmt)$/ || $bare=~/^(?:from|to)_/ || $cat=~m{^(?:conversion|io/format)$}) ? 1 : 0 }
sub m32row { my ($a,$p)=@_; my %b=%$a; my $zb=$b{crate} eq 'zb'; my $ia=$zb?'Err(BoundError::InvalidArgument)':'Error::InvalidArgument';
  my $MX='MAX-timeout: a run past the 10-minute budget is recorded SKIPPED_AFTER_TIMEOUT, which counts as FAIL in the verdict space and in pm-coverage and is never omitted';
  my $red="Red until $ZFX{'ZF-M32'} lands and Rusty resolves the pin ($PMNOTE): the pinned source predates M32 (MAX_PRECISION there is the pre-M32 value in the CSV header) and no byte budget exists (common/buf.rs:55-70 calls try_reserve_exact only).";
  if($p eq 'MAX+1'){ $b{class}='prec-PM_MAX_P+1'; $b{class_def}="p = PM_MAX_P + 1, one above the named pin ($PMNOTE)";
    $b{required_result}="refusal: $ia raised up front by the precision check (zf p_assertion num.rs:55-57; zb precision.rs:17), before any allocation; any returned value is FAIL; a timeout is FAIL".($zb?" (passed as usize, precision.rs:9; the CONTRACT u32 text is fixed by $ZFX{'ZB-D05'})":'').". $red $MX";
    $b{refusal_variant}=$ia; $b{oracle_x86_64}="spec: refusal $ia (M32-DECISION.md); no numeric gold"; return \%b }
  $b{class}='prec-PM_MAX_P'; $b{class_def}="p = PM_MAX_P, the named MAX_PRECISION pin ($PMNOTE); a single buffer above 2^31 B refuses, exactly 2^31 B is accepted";
  if(m32basic($b{family},$b{category})){
    $b{required_result}="correct result at p = PM_MAX_P in ToEven".($zb?' (enclosure containing the exact value, as in the normal row)':' (as in the normal row)').", bit-identical to the one gold on every active target (P1). Basic op: Paul M32 ruling 2026-10-30 7:18 PM CT requires add, sub, mul, div, sqrt and every conversion to be correct at MAX in zf and zb, so a refusal is FAIL; a crash is FAIL; an OOM kill is FAIL; a timeout is FAIL. $red $MX";
    $b{refusal_variant}='' }
  else { $b{required_result}="TRACE-PENDING [RC-FAIL]: no outcome is pinned at p = PM_MAX_P until the cross-file allocation trace (design s28.4) compares the largest reachable single buffer against 2^31 B; every run of this row is FAIL, never PASS, until the trace pins it. $red $MX";
    $b{refusal_variant}=''; $b{oracle_x86_64}='TRACE-PENDING: no gold until the design s28.4 trace pins the outcome' }
  \%b }
sub row { my %a=@_; pin1(\%a); pin2(\%a); pin3(\%a); my @R=(\%a); for my $t (@TWIN){ $$t{oracle_i686}=$a{oracle_i686}; push @R,$t } @TWIN=();
  if($a{crate} eq 'zf' && $a{class} eq 'normal' && ($a{verdict_space}//'') ne 'NA_PROOF' && $a{family}=~$M30RE && !$M30DONE{$a{family}}++){ for my $e (24,27){ push @R,{%a, class=>"m30-p2^$e", class_def=>"p = 2^$e on both targets (ZF-M30 cost-overflow threshold, an input point; P1 one gold)", inputs=>'the normal operand of the family', modes=>'ToEven', tier=>'heavy only (own shard pm-max, run on both targets; budget 10 min per run; timeout = FAIL)', required_result=>"correct result at p = 2^$e in ToEven, bit-identical to the one gold on both targets (P1). Each of these is FAIL, in debug and release builds: a panic, a wrapped cost estimate, an overflowed cost estimate, a refusal, a capped precision; MAX-timeout skip applies (design s4). Red until $M30 lands: at 7955836 the series and reduction cost estimates overflow at this p on some targets (series.rs:156,238,244; tan.rs:128; common/util.rs:55-64 calc_mul_cost; cos.rs:35 reduction_cost).", refusal_variant=>'', verdict_space=>'PASS|FAIL|SKIPPED_AFTER_TIMEOUT'} } }
  @R=grep { pin4($_) } @R; @R=map { $$_{class} eq 'prec-over-MAX' ? ($_, m32row($_,'MAX'), m32row($_,'MAX+1')) : $_ } @R;
  for my $r (@R){ my %a=%$r; { my @t=pm_target_specific(@a{@F17}); push @LEFT,"target(@t) $a{family} $a{class}" if @t; } { my @h=pm_hedges($a{required_result},$a{refusal_variant}); push @LEFT,"hedge(@h) $a{family} $a{class}: $a{required_result} || $a{refusal_variant}" if @h; } push @out,[map { defined $_ ? s/\boracle_x86_64\b/the oracle column/gr : $_ } @a{@F17}]; $cnt{"$a{crate} $a{verdict_space}"}++ } }
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
  my $V = sub { my $v=shift; return "Err(BoundError::$v)" if $crate eq 'zb'; return 'None' if $ret eq 'option'; return "Err(Error::$v)" if $ret eq 'result'; return "false (documented refusal, $v-class input)" if $ret eq 'bool'; "NaN whose error kind is Error::$v" };
  my $proofp="docs/test-plan/pm-proofs/$crate/".($key=~s/[^A-Za-z0-9_:.-]/_/gr).'.md';
  my $gold = sub { my $g=($goldfor{$crate}{lc $bare}||[]); @$g ? join(' | ',map {"gold:$_"} @$g) : "gold:".($crate eq 'zb'?'zenith-bound/tests/goldens/pm':'zenith-float-num/tests/pm/golds')."/$crate-".($key=~s/[^A-Za-z0-9_]+/_/gr).".jsonl\@sha256:UNGENERATED (pm-coverage FAIL until written by PM-3 and pinned in SHA256SUMS)" };
  my $ox = $crate eq 'zb' ? "primary: mpmath 1.3.0 at prec+256 via goldlib (zenith-bound/tools, ORACLE_PY=$ENV); second check: python-flint 0.9.0 arb, ".$gold->()
        : ($mpfrnative{$bare} && $key=~/^ExactNum::/) ? "primary: $MPFR"
        : ($grp eq 'G13 random (feature=random, R2)') ? "primary: mpmath 1.3.0 in $ENV generates the gold (the oracle script replays the splitmix64 seed stream and builds each value at p bits); the row checks a digest of it: digest:".$gold->()
        : ($bare=~/gcd|lcm|mod_|egcd|pollard|miller|factor|binomial|isqrt|int/i) ? "primary: mpmath 1.3.0 exact integer arithmetic (gmpy2 backend) in $ENV; second check: python-flint 0.9.0, ".$gold->()
        : "primary: mpmath 1.3.0 Ziv pair p+256/p+512 must agree, in $ENV, ".$gold->();
  my $oi = $crate eq 'zb' ? 'same pinned gold file + SHA as x86_64; byte-identical BallRow digest to the x86_64 run'
        : 'gold:zenith-float-num/tests/pm/golds/i686/'.($key=~s/[^A-Za-z0-9_]+/_/gr).".jsonl\@sha256:UNGENERATED (cut at PM-3 from the x86_64 run checked against the same primary oracle as oracle_x86_64; pm-coverage FAIL until pinned)";
  my $chk = $crate eq 'zb' ? 'containment: exact value in [mid-rad, mid+rad]; rad <= tightness baseline' : 'exact limbs + sign + inexact flag + error kind equal the pinned oracle';
  my %base=(crate=>$crate,family=>$key,source=>$src,category=>$cat,group=>$grp,oracle_x86_64=>$ox,oracle_i686=>$oi,check=>$chk,proof=>'',refusal_variant=>'',_ret=>$ret,_hasrm=>$hasrm);
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
    else { my $t=$TIE{$key} or die "TIE RESOLUTION MISSING $key (inputs/tie-resolution.tsv)\n"; my ($res,$tm,$basis,$note)=@$t; my $TP='constructed ties at p in {32,64,128}';
      my %DIR=(ToEven=>'the even neighbour',ToZero=>'the neighbour nearer zero',FromZero=>'the neighbour farther from zero');
      if($res eq 'r4'){ $R->('tie','exact half-ulp ties at p (fixed internal mode; outcome from the R4 row)',$TP,'the mode the R4 row pins','fast+heavy','','',r4class=>'tie') }
      elsif($res eq 'exact'){ $R->('tie','exact half-ulp ties at p (no rounding step)',$TP,'none (no rounding step)','fast+heavy',"exact result: the family has no rounding step (source $src; $basis), so the constructed tie input yields its exact value, bit-identical to the one gold (P1)",'',keep=>1) }
      elsif($res eq 'ieee'){ $R->('tie','exact half-ulp ties at the IEEE format precision (fixed roundTiesToEven)','constructed ties at the format precision (binary32: 24 bits, binary64: 53 bits)','ToEven (IEEE roundTiesToEven, fixed by the format)','fast+heavy',"correctly rounded in IEEE 754-2019 roundTiesToEven at the format precision (binary32: 24 bits, binary64: 53 bits): the tie goes to the even significand, bit-exact against the IEEE oracle ($basis; source $src). The format fixes the precision, so the p axis and ZF-M27 do not apply",'',keep=>1) }
      elsif($res eq 'mode'){ $R->('tie',"exact half-ulp ties at p (fixed internal mode $tm)",$TP,"$tm at p in {32,64,128} (fixed internal mode, $src)",'fast+heavy',"correctly rounded at p in $tm, the mode the family fixes internally ($src): the tie goes to $DIR{$tm}, rounded once at rp(p,64), bit-identical to the one gold (R3, $M27). Red until ZF-M27 lands: the p=32 cell is FAIL until rounding happens at rp(p,64).",'',keep=>1) }
      elsif($res eq 'none'){ $R->('tie','exact half-ulp ties at p (fixed internal mode None)',$TP,"None at p in {32,64,128} (fixed internal mode, $src)",'fast+heavy',"R6: None rounds as RNDZ at rp(p,64), so the tie goes to the neighbour nearer zero at rp(p,64) bits, bit-identical to the one gold (R6, P1; source $src)",'',keep=>1) }
      elsif($res eq 'rm'){ $R->('tie','exact half-ulp ties at p (family takes a RoundingMode: '.(index($key,"Context::")==0?'ctx.rm':'rm parameter of the generated wrapper').", $src)",'constructed ties',$X7,$FX,'ToEven -> even neighbour, ToOdd -> odd neighbour, directed modes per direction, None -> RNDZ@rp(p,64) (R6)','',keep=>1) }
      elsif($res eq 'pinned-red'){ $R->('tie',"exact half-ulp ties at p (outcome pinned: $tm; $basis)",$TP,"$tm (pinned, $basis)",'fast+heavy',"$note, bit-identical to the one gold (P1). [RUSTY-CONFIRM] Red until Rusty confirms this against the source $src (owner Rusty; $basis): every run is FAIL until then",'',keep=>1) }
      elsif($res eq 'red'){ $R->('tie','exact half-ulp ties at p (outcome unpinned)',$TP,'unpinned (fixed internal mode not stated)','fast+heavy',"[PM1-UNPINNED] FAIL until Rusty states $note. The single outcome is pinned when that fact lands (source $src)",'',keep=>1) }
      else { die "TIE RESOLUTION BAD $res for $key\n" } } }
  # 5: exact range thresholds (EXPONENT_MIN/MAX, not i64)
  my $ovf = $crate eq 'zb' ? 'enclosure, or the all-reals ball (rad=Infinity) on midpoint overflow (CONTRACT 10.4); never a finite ball excluding the true value' : 'ExponentOverflow(sign): ±Inf for modes away from zero; max finite flagged inexact for modes toward zero (ToZero; Down for +; Up for -) per CONTRACT 10.4 description of zf';
  my $unf = $crate eq 'zb' ? 'radius absorbs the lost value and is never zero (CONTRACT 10.4, 3.2)' : 'correctly rounded subnormal / signed zero per mode (EXPONENT_MIN gradual underflow)';
  $R->('ovf-max-finite','result = max finite at p (exponent EXPONENT_MAX = 2^62-1)','operands whose exact result is max finite',$hasrm?$X7:$M7,$FX,$crate eq 'zb'?$ovf:'exact max finite','',r4class=>'huge-exp');
  $R->('ovf-max-plus-half-ulp','exact result = max finite + half ulp (and just below/above)','constructed',$hasrm?$X7:$M7,$FX,$ovf,$crate eq 'zb'?'Ok(all-reals ball)':'ExponentOverflow(Sign)',r4class=>'huge-exp');
  $R->('unf-min-normal','result = smallest normal 2^(EXPONENT_MIN-1)','constructed',$hasrm?$X7:$M7,$FX,$unf,'',r4class=>'tiny-exp');
  $R->('unf-min-subnormal','result = smallest subnormal 2^(EXPONENT_MIN-p)','constructed',$hasrm?$X7:$M7,$FX,$unf,'',r4class=>'subnormal-like');
  $R->('unf-half-min-subnormal-tie','exact result = half the smallest subnormal (tie to 0 or min subnormal)','constructed',$hasrm?$X7:$M7,$FX,$crate eq 'zb'?$unf:'ToEven, ToZero, Down(+) -> +0; Up(+), FromZero -> min subnormal; ToOdd -> min subnormal; None -> RNDZ (+0)','',r4class=>'tiny-exp');
  if($i64x{$bare}){ $R->('i64-extreme-exponent-arg','exponent argument i64::MIN, i64::MIN+1, i64::MAX-1, i64::MAX, EXPONENT_MAX±1, EXPONENT_MIN±1','constructed',$hasrm?$X7:$M7,$FX,'no wrap; an exponent above EXPONENT_MAX gives ExponentOverflow(sign), one below EXPONENT_MIN gives the rounded underflow value in rm; ilogb and frexp are exact','ExponentOverflow(Sign)',keep=>1) }
  # 4: huge trig args
  if($trig{$bare}){ my $hz=$hasrm?$X7:$M7; my $FM='(FIXLIST:15 ZF-01, FIXLIST:127 ZF-M22)';
    for my $e ('2^16-1','2^16','2^16+1'){ my $what = $bare eq 'rem_pi' ? 'x mod pi correctly rounded at p in rm against the mpmath value (FIXLIST:31 ZF-17)' : ($crate eq 'zb'?'an enclosure of the mpmath value':'the correctly rounded mpmath value at p in rm');
      $R->("trig-huge-arg e=$e","|x| with exponent $e (a ZF-M22 regression point; no refusal: the result is required)",'x = (1+2^-k)*2^e, random mantissas',$hz,'fast(ToEven,p=64)+heavy',"$what; identical on i686 and x86_64; a refusal is FAIL, a wrong value is FAIL (Spence pmrev5 timing: up to 0.07 s per row at p=1025, about 0.002 s at e=2^16+-1) $FM",'',keep=>1) }
    my $MA=$V->('MemoryAllocation');
    $R->('trig-huge-arg e=EXPONENT_MAX','|x| with exponent EXPONENT_MAX = 2^62-1 (defs.rs:49)','x = (1+2^-k)*2^EXPONENT_MAX',$hz,'fast(ToEven,p=64)+heavy',"refusal: $MA, raised before allocating, because the reduction needs more pi bits than $PIB; identical on i686 and x86_64; any returned value is FAIL $FM",$MA,keep=>1);
    $R->('trig-pi-budget-under','|x| whose reduction needs pi bits just under the pi-bit budget constant','constructed from the published constant',$hz,'heavy',$UNP."the correctly rounded mpmath value at p in rm; FAIL until $PIB is published and the gold is cut $FM",'',keep=>1);
    $R->('trig-pi-budget-over','|x| whose reduction needs pi bits just over the pi-bit budget constant','constructed from the published constant',$hz,'heavy',$UNP."refusal $MA raised before allocating; FAIL until $PIB is published and the gold is cut $FM",$MA,keep=>1) }
  # 12: T4 mul operand pairs
  if($bare eq 'mul' && $crate eq 'zf'){ for my $w (69,70,71){ my $l=2*$w; $R->("mul-pair-T4 sm=$w lg=$l","unbalanced operands: sm=$w words, lg=$l words (mantissa/mul.rs:53 lg.len()/2 >= sm.len() && sm.len() > 70)","p_a=$w*W, p_b=$l*W, both orders",$hasrm?$X7:$M7,$FX,'correctly rounded at p in rm, equal across the switch','',keep=>1) } }
  # primality (Paul)
  if($crate eq 'zf' && $key=~/^(?:ball::ziv_round_vec|ExactNum::meijer_g|ops::special::beta)$/){ $R->('doc-intra-links','rustdoc intra-doc links in the docs of this family resolve ('.$ZFX{'ZF-D08'}.')','docs at '.$src,'n/a','fast',"FAIL until $ZFX{'ZF-D08'} lands: cargo doc with RUSTDOCFLAGS=\"-D warnings\" must report no broken intra-doc link in the docs of $bare; FIXLIST ZF-D08 names it as broken at 7955836",'',keep=>1) }
  if($bare=~/prime|primality|miller/ && $bare!~/^(ai|bi)_prime$/){ $R->('doc-probable-prime','primality test naming/doc rule (Paul)','docs at '.$src,'n/a','fast','name or docs say "probable prime", or the function is deterministic for all inputs; currently miller_rabin docs (modular.rs:105-108) say neither -> FAIL until FIXLIST doc fix','',keep=>1) }
  # pmrev4 item 4: rows that carry the single pinned refusal of the families whose value rows now forbid it
  if($crate eq 'zf' && $key=~/^roots::(brent|illinois)$/){ $R->('no-sign-change','bracket [a,b] with f(a) and f(b) of the same sign (no sign change)','constructed brackets',$M7,$FX,'refusal: NaN whose error kind is Error::InvalidArgument, held in the returned Option (roots.rs:238, 407 return Option<ExactNum>); None is FAIL; a returned root is FAIL','NaN whose error kind is Error::InvalidArgument (held in the returned Option)',keep=>1) }
  if($crate eq 'zf' && $key eq 'ode::rk4'){ $R->('bad-step-interval','step h <= 0, non-finite h, zero step count, empty interval, reversed interval','constructed',$M7,$FX,'refusal: Error::InvalidArgument; None is FAIL; any returned value is FAIL; rk4 returns Option at 7955836 (ode.rs:73), so this row FAILs until rk4 reports the kind: FIXLIST API row needed','Error::InvalidArgument',keep=>1) }
  if($crate eq 'zf' && $key=~/^ExactNumArray::(lu|qr|svd|eigen)_decomp$/){ my %dg=(lu=>'an exactly singular matrix (a zero pivot after partial pivoting; exact ties go to the lowest index)',qr=>'a zero R_ii (rank-deficient column)',svd=>'repeated singular values',eigen=>'repeated eigenvalues'); my $d=$dg{$1};
    $R->('degenerate',"matrix with $d",'constructed exact matrices',$M7,$FX,"None (pinned per row in the gold for $d); Some(factors) is FAIL",'None',keep=>1) }
  if($crate eq 'zf' && $key eq 'ExactNum::from_raw_parts'){ $R->('raw-parts-invalid','parts invalid at p: exponent outside [EXPONENT_MIN, EXPONENT_MAX], limb count not matching p, p outside [1, MAX_PRECISION]','constructed raw parts',$M7,$FX,'refusal: '.$V->('InvalidArgument').'; any returned value is FAIL',$V->('InvalidArgument'),keep=>1) }
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

{ my $DER="Derivation: every zf heap buffer is a Vec<Word> with Word = u64 (defs.rs:16; WordBuf::new common/buf.rs:55-70 reserves sz words), so a buffer is 8*sz bytes: 2^31-1 and 2^31+1 are not multiples of 8 and no input produces them; 2^31 = WordBuf::new(2^28) is reachable only crate-internally. No public precision reaches 2^31 B exactly: a mantissa at p = 2^32-1 (the M32-DECISION.md value) is at most 2^26 words = 2^29 B (base.rs:51-53), a product at most 2^27 words = 2^30 B (base.rs:502), and the fft_mul buffer jumps from 1612234800 B (L1+L2 = 33546240 words) to 3222896688 B (33546241 words) (fft.rs:437-443, fft_params fft.rs:209-222, sqrt_int common/util.rs:42-52; derivation in design s28.3). So these rows call the M32 byte check with the requested size n in bytes. 'Larger than 2^31 bytes' (M32-DECISION.md) is read literally: 2^31 is not larger, so it is accepted.";
  for my $bb (['2^31-1','Ok: accepted (n is not larger than 2^31)',''],['2^31','Ok: accepted (n = 2^31 is not larger than 2^31, read literally)',''],['2^31+1','refusal: Error::MemoryAllocation from the check, before try_reserve_exact; Ok is FAIL','Error::MemoryAllocation']){
    my %a=(crate=>'zf',family=>'M32 byte budget check (WordBuf::new)',source=>'zenith-float-num/src/common/buf.rs:55',category=>'limits',group=>'M32-bytes',class=>"m32-bytes-$$bb[0]",class_def=>"a single buffer request of n = $$bb[0] bytes (M32-DECISION.md byte budget)",inputs=>"n = $$bb[0] bytes passed to the M32 up-front size check. $DER",modes=>'n/a (no rounding)',tier=>'fast (unit test of the check)',required_result=>"$$bb[1]. Red until $ZFX{'ZF-M32'} lands: at 7955836 there is no byte check (buf.rs:55-70 calls try_reserve_exact only).",refusal_variant=>$$bb[2],oracle_x86_64=>'spec: M32-DECISION.md (sha256 e6bcf204) byte budget; no numeric gold',oracle_i686=>"one result for every active target in docs/test-plan/inputs/targets.csv (P1); a pending target is reported NOT-RUN, never PASS",check=>'unit test: assert the check result',verdict_space=>'PASS|FAIL',proof=>'');
    push @out,[@a{@F17}]; $cnt{"zf PASS|FAIL"}++ } }
die "pin1: ".scalar(@LEFT)." rows left with either/or text:\n".join("\n",@LEFT[0..($#LEFT<40?$#LEFT:39)])."\n" if @LEFT;
open my $o,'>',$OUT or die; my @H=qw(crate family source category group class class_def inputs modes tier required_result refusal_variant oracle gold_scope check verdict_space proof);
sub csvq { my $s=shift; $s="" unless defined $s; $s=~s/\s*\n\s*/ /g; $s=~s/"/""/g; "\"$s\"" }
print $o "# precision-matrix-edges v12 (Harper 2026-10-30 pmrev12: PM_MAX_P pin, Paul M32 ruling 7:18 PM CT; P1 bit-identical everywhere and P6 64-bit minimum, pop-logs/PRINCIPLES.md; ZF-M27..M32; targets in docs/test-plan/inputs/targets.csv, no per-target columns). Lines starting with # are metadata; the CSV header is the first line without #.\n# source commit $SRC_COMMIT (generator refuses if zenith-float-num/{src,tests} or zenith-bound/{src,tests} differ, except files sha256-pinned in gen_edges.pl %SRCPIN)\n";
print $o "# input docs/test-plan/inputs/$_ sha256:$PIN{$_}\n" for sort keys %PIN;
print $o "# PM_MAX_P = zf MAX_PRECISION = $$PMX{value} ($$PMX{form}) from $$PMX{file}:$$PMX{line} at $SRC_COMMIT (file sha256:$$PMX{sha}); M32 in source: ".($$PMX{m32}?"yes":"no")."; expected $$PMX{expected} at $$PMX{expected_commit} $$PMX{expected_at}; ".($$PMX{m32}&&$$PMX{match}?"resolved":"OPEN: RC FAILs until the pinned source includes M32 and matches")."\n";
print $o join(',',map {csvq($_)} @H),"\n"; print $o join(',',map {csvq($_)} @$_),"\n" for @out; close $o;
my %fc; $fc{$$_[0]}++ for @fam; my %wp=(); $wp{$$_[1]}=1 for grep {$$_[6]=~/working-precision/} @fam;
print "families zf=$fc{zf} zb=$fc{zb} rows=",scalar(@out),"\n"; print "$_ $cnt{$_}\n" for sort keys %cnt; print "wp families: ",join(' ',sort keys %wp),"\n";
my %oth; $oth{$$_[2]}++ for grep {$$_[0] eq 'zf' && $$_[5] eq 'other'} @fam; print "other keys ",scalar(keys %oth)," bare; r4 fams ",scalar(keys %r4fam),"; missing: ",join(' ',grep {!$r4fam{$_}} sort keys %oth),"\n";

{ my %sh; for my $r (@out){ my ($cr,$fm,$tier,$vs)=@{$r}[0,1,9,15]; my $k = $vs eq 'NA_PROOF' ? 'proof' : $tier=~/pm-max-i686/ ? 'pm-max-i686' : $tier=~/pm-max-x86_64/ ? 'pm-max-x86_64' : $tier=~/shard pm-max\b/ ? 'pm-max' : $tier=~/fast/ ? 'fast' : 'heavy'; $sh{"$k|$cr|$fm"}++ }
  open my $so,'>',$SHOUT or die; print $so "shard,crate,family,rows\n"; for (sort keys %sh){ my ($a,$b,$c)=split /\|/; print $so join(',',map {csqq($_)} $a,$b,$c,$sh{$_}),"\n" } close $so; print "shard lines: ",scalar(keys %sh),"\n" }
sub csqq { my $s=shift; $s=~s/"/""/g; "\"$s\"" }
open my $n5,'>',$N5OUT or die; print $n5 "crate,family,N,citations\n"; my %fam5;
for my $F (@fam){ my ($crate,$key)=@$F; my $c=$n5{$key}; print $n5 join(',',map {csvq($_)} $crate,$key,$c?5:3,$c?join(' | ',@$c):'none (N=3)'),"\n"; $fam5{$c?5:3}++ } close $n5;
print "N5 families: $fam5{5}; N3 families: $fam5{3}\n";
for my $i (0..$#N5){ my @miss=grep {$_ ne '*' && !$n5hit{$i}{$_}} @{$N5[$i][2]}; print "N5 unmatched names (entry $i): @miss\n" if @miss }
for my $p (@PROOFS){ my ($pp,$key,$src,$bare)=@$p; my $f="$T/$pp"; next if -e $f; open my $pf,'>',$f or die "$f"; print $pf "# NA_PROOF: $key tie row\n\nClaim: `$key` ($src) takes no RoundingMode and performs no rounding (Spence G1-G3: accessor, construct-copy or exact integer/modular), so an exact half-ulp tie at p cannot arise.\n\nEvidence: TO BE WRITTEN by the proving agent (cite the signature and body at $src).\n\nSecond-agent sign-off: UNSIGNED\n\n<!-- pm-coverage FAILs this row until the line above reads `Second-agent sign-off: <agent> <date> <commit>` -->\n"; close $pf }
print "dropped ($_): $DROP{$_}\n" for sort keys %DROP;
print "proof placeholders: ",scalar(@PROOFS),"\n";

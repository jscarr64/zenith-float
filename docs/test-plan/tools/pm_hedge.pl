# Harper 2026-10-30 (Spence pmrev5 item 1): the one hedge rule, shared by tools/gen_edges.pl (refuses to write) and tools/lint.pl (BAD).
# Applies to required_result (col 10) and refusal_variant (col 11) of every row. pm_hedges($c10,$c11) returns the broken rule names.
#  UNCLEAR     the word UNCLEAR
#  TBD/TODO    \bTBD\b, \bTODO\b (any case)
#  either/or   \bor\b, \beither\b (any case), a standalone '|'
#  else        \belse\b (any case; covers ', else X')
#  unless      \bunless\b (any case)
#  otherwise   \botherwise\b (any case)
#  slash       a '/' between alternatives. Each whitespace token holding '/' is stripped of leading ([{"'` and trailing
#              )]}"'`,;:. and is allowed only if it is (P) a file path: segments of [~\w.-] joined by '/', ending in
#              .rs .md .csv .tsv .py .txt .toml .jsonl .log .pl with an optional :lines suffix, or (M) math: every '/'
#              in it has a digit right before or right after it, or the token has one of = * ^ in it, or (X) every '/' in the token sits inside
#              parentheses (an expression such as cos(pi/n)). Every other '/' (a standalone '/', word/word, Mode/Mode,
#              Down(+)/FromZero) is a hedge.
#  rv-multi    refusal_variant holding ';', or more than one distinct outcome token
#              (Error::X, BoundError::X, None, Some, Ok) - one row, one refusal.
sub pm_slash_ok { my $t=shift; $t=~s/^[(\[{"'`]+//; $t=~s/[)\]}"'`,;:.]+$//;
  return 1 if $t=~m{^[~\w.\-]+(?:/[~\w.\-]+)+\.(?:rs|md|csv|tsv|py|txt|toml|jsonl|log|pl)(?::[\d,\-]+)?$};
  return 1 if $t=~/[=*^]/;
  { my ($d,$in)=(0,1); for my $ch (split //,$t){ $d++ if $ch eq '('; $d-- if $ch eq ')'; $in=0 if $ch eq '/' && $d<=0 } return 1 if $in }
  my @p=split m{/},$t,-1; return 0 if @p<2;
  for my $i (0..$#p-1){ return 0 unless $p[$i]=~/\d$/ || $p[$i+1]=~/^\d/ } 1 }
sub pm_hedges { my ($a,$b)=@_; $a//=''; $b//=''; my $t="$a $b"; my @r;
  push @r,'UNCLEAR' if $t=~/UNCLEAR/;
  push @r,'TBD/TODO' if $t=~/\b(?:TBD|TODO)\b/i;
  push @r,'either/or' if $t=~/\b(?:or|either)\b/i || $t=~/(?:^|\s)\|(?:\s|$)/;
  push @r,'else' if $t=~/\belse\b/i;
  push @r,'unless' if $t=~/\bunless\b/i;
  push @r,'otherwise' if $t=~/\botherwise\b/i;
  push @r,'slash' if grep { m{/} && !pm_slash_ok($_) } split /\s+/,$t;
  my %o; $o{$1}=1 while $b=~/((?:Bound)?Error::\w+|\bNone\b|\bSome\b|\bOk\b)/g;
  push @r,'rv-multi' if $b=~/;/ || keys(%o)>1;
  @r }
1;

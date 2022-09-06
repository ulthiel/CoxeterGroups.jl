// This script benchmarks the speed at which Magma can multiply Coxeter group elements.
// The script will compute all elements of the given group, then multiply all ordered pairs of elements in the group.
// In order to try to account for the interpretation overhead in Magma, we run another loop over all ordered pairs which
// does a very simple computation, and remove that from the time.

SetColumns(0);

function AllGroupEelements(W)
    queue := {@ W.0 @};
    gens := Generators(W);
    pos := 1;
    while pos le #queue do
        w := queue[pos];
        pos +:= 1;

        for g in gens do
            Include(~queue, w * g);
            Include(~queue, w * g^-1);
        end for;
    end while;

    return queue;
end function;

function SumLengths(W, Welts)
    one := W.0;
    sum := 0;
    for x in Welts, y in Welts do
        sum +:= x eq one select 1 else 0;
    end for;
    return sum;
end function;

function Multiply(W, Welts)
    one := W.0;
    sum := 0;
    for x in Welts, y in Welts do
        sum +:= (x * y) eq one select 1 else 0;
    end for;
    return sum;
end function;


procedure PrintMultiplyStats(benchName, W)
    Welts := AllGroupEelements(W);

    start := Realtime();
    _ := SumLengths(W, Welts);
    sumtime := Realtime(start);

    start := Realtime();
    _ := Multiply(W, Welts);
    multtime := Realtime(start);

    usPerElt := (multtime - sumtime) / #Welts^2 * 1000 * 1000;

    printf "%o: All-pairs multiplication (%o elements), about %o μs per element.\n",
        benchName, #Welts, usPerElt;
end procedure;

PrintMultiplyStats("GrpFPCox A5", CoxeterGroup(GrpFPCox, "A5"));
PrintMultiplyStats("GrpPermCox A5", CoxeterGroup(GrpPermCox, "A5"));
PrintMultiplyStats("SymmetricGroup A5", SymmetricGroup(6));

print "";
PrintMultiplyStats("GrpFPCox A6", CoxeterGroup(GrpFPCox, "A6"));
PrintMultiplyStats("GrpPermCox A6", CoxeterGroup(GrpPermCox, "A6"));
PrintMultiplyStats("SymmetricGroup A6", SymmetricGroup(7));


quit;

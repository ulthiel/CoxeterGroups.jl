// The purpose of this script is to try to get an idea of the speed of Magma's GrpFPCox implementation.
// It should be run like this:
//     magma -b type:=A6 magma.m
//
// The script will compute all elements of the given group, then multiply all ordered pairs of elements in the group.
// In order to try to account for the interpretation overhead in Magma, we run another loop over all ordered pairs which
// does a very simple computation, and remove that from the time.

SetColumns(0);

function AllElementsToLength(W : maxLength := -1)
    queue := {@ W.0 @};
    pos := 1;
    while pos le #queue do
        w := queue[pos];
        pos +:= 1;

        for s := 1 to Rank(W) do
            ws := w * W.s;
            if #w lt #ws and (maxLength eq -1 or #ws le maxLength) then
                Include(~queue, ws);
            end if;
        end for;
    end while;

    return queue;
end function;

function SumLengths(Welts)
    sum := 0;
    for x in Welts, y in Welts do
        sum +:= #x + #y;
    end for;
    return sum;
end function;

function Multiply(Welts)
    sum := 0;
    for x in Welts, y in Welts do
        sum +:= #(x * y);
    end for;
    return sum;
end function;

W := CoxeterGroup(GrpFPCox, type);
Welts := AllElementsToLength(W);

printf "Working in type %o with %o elements\n", type, #Welts;

start := Realtime();
SumLengths(Welts);
sumtime := Realtime(start);

start := Realtime();
Multiply(Welts);
multtime := Realtime(start);

usPerElt := (multtime - sumtime) / #Welts^2 * 1000 * 1000;

printf "It took %o seconds to sum the lengths of all pairs of elements,\n", sumtime;
printf "and     %o seconds to sum the lengths of the product of all pairs of elements.\n", multtime;
printf "Hence the internal Magma multiplication takes about %o μs per element.\n", usPerElt;




quit;

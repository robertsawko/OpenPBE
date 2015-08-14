# header
BEGIN {
    Iteration=0
    resetCounters()
}

# reset counters used for variable postfix
function resetCounters() {
    CourantMaxCnt=0
    CourantMeanCnt=0
    e_oilCnt=0
    e_oilFinalResCnt=0
    e_oilItersCnt=0
    epsilonmCnt=0
    epsilonmFinalResCnt=0
    epsilonmItersCnt=0
    e_waterCnt=0
    e_waterFinalResCnt=0
    e_waterItersCnt=0
    executionTimeCnt=0
    kmCnt=0
    kmFinalResCnt=0
    kmItersCnt=0
    pCnt=0
    pFinalResCnt=0
    pItersCnt=0
    SeparatorCnt=0
    TimeCnt=0
    # Reset counters for general Solving for extraction
    for (varName in subIter)
    {
        subIter[varName]=0
    }
}

# Extract value after columnSel
function extract(inLine,columnSel,outVar,a,b)
{
    a=index(inLine, columnSel)
    b=length(columnSel)
    split(substr(inLine, a+b),outVar)
    gsub("[,:]","",outVar[1])
}

# Iteration separator (increments 'Iteration')
/^[ \t]*Time = / {
    Iteration++
    resetCounters()
}

# Time extraction (sets 'Time')
/^[ \t]*Time = / {
    extract($0, "Time = ", val)
    Time=val[1]
}

# Skip whole line with singularity variable
/solution singularity/ {
    next;
}

    # Extraction of any solved for variable
    /Solving for/ {
        extract($0, "Solving for ", varNameVal)

        varName=varNameVal[1]
        file=varName "_" subIter[varName]++
        file="./logs/" file
        extract($0, "Initial residual = ", val)
        print Time "\t" val[1] > file

        varName=varNameVal[1] "FinalRes"
        file=varName "_" subIter[varName]++
        file="./logs/" file
        extract($0, "Final residual = ", val)
        print Time "\t" val[1] > file

        varName=varNameVal[1] "Iters"
        file=varName "_" subIter[varName]++
        file="./logs/" file
        extract($0, "No Iterations ", val)
        print Time "\t" val[1] > file
    }

# Extraction of CourantMax
/Courant Number / {
    extract($0, "max: ", val)
    file="./logs/CourantMax_" CourantMaxCnt
    print Time "\t" val[1] > file
    CourantMaxCnt++
}

# Extraction of CourantMean
/Courant Number / {
    extract($0, "mean: ", val)
    file="./logs/CourantMean_" CourantMeanCnt
    print Time "\t" val[1] > file
    CourantMeanCnt++
}

# Extraction of executionTime
/ExecutionTime = / {
    extract($0, "ExecutionTime = ", val)
    file="./logs/executionTime_" executionTimeCnt
    print Time "\t" val[1] > file
    executionTimeCnt++
}

# Extraction of Separator
/^[ \t]*Time = / {
    extract($0, "Time = ", val)
    file="./logs/Separator_" SeparatorCnt
    print Time "\t" val[1] > file
    SeparatorCnt++
}

# Extraction of Time
/^[ \t]*Time = / {
    extract($0, "Time = ", val)
    file="./logs/Time_" TimeCnt
    print Time "\t" val[1] > file
    TimeCnt++
}


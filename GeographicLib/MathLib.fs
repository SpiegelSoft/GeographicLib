namespace GeographicLib

open System

module MathLib =
    let atanh(x) = (log(1.0 + x) - log(1.0 - x))/2.0
    let rec getMachineEpsilon eps = if 1.0 + eps = 1.0 then eps else getMachineEpsilon ((abs eps) / 2.0)
    let eatanhe(x, es) = if es >= 0.0 then es * atanh(es * x) else -es * atan(es * x)
    let sq(x : float<'a>) = x * x
    let minFloat = Double.Epsilon
    let epsilon = getMachineEpsilon 0.1
    let epsilonSquared = sq(epsilon)
    let epsilonToTheFourth = sq(epsilonSquared)
    let tolerance = sqrt epsilon
    let tolerance0 = tolerance * epsilon |> sqrt |> sqrt
    let piOver2 = Math.PI / 2.0
    let polyval(N, p : float[], x) =
        let mutable y = if N < 0 then 0.0 else p.[0]
        let mutable N = N
        let mutable i = 1
        while N > 0 do
            N <- N - 1
            y <- y * x + p.[i]
            i <- i + 1
        y


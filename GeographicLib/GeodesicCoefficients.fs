namespace GeographicLib

open System

module GeodesicCoefficients =
    let GeodesicOrder =
        match Constants.Precision with
        | 2 -> 6
        | 1 -> 3
        | 3 -> 7
        | _ -> 8
    let nA1 = GeodesicOrder
    let nC1 = GeodesicOrder
    let nC1p = GeodesicOrder
    let nA2 = GeodesicOrder
    let nC2 = GeodesicOrder
    let nA3 = GeodesicOrder
    let nA3x = nA3
    let nC3 = GeodesicOrder
    let nC3x = (nC3 * (nC3 - 1)) / 2;
    let nC4 = GeodesicOrder
    let nC4x = (nC4 * (nC4 + 1)) / 2;
    // Size for temporary array
    // nC = max(max(nC1_, nC1p_, nC2_) + 1, max(nC3_, nC4_))
    let nC = GeodesicOrder + 1
    [<Literal>]
    let maxit1 = 20
    
    let Fourier n (coeff : float[]) eps =
        let eps2 = sqrt eps
        let mutable d = eps
        let mutable o = 0
        [|1..n|] |> Array.map (fun l ->
            let m = (n - l) / 2
            let v = d * MathLib.polyval(m, coeff.[o..], eps2) / coeff.[o + m + 1]
            o <- o + m + 2
            d <- d * eps
            v) |> Array.append [|0.0|]

    // The coefficients C1[l] in the Fourier expansion of B1
    let C1Fourier eps =
        let coeff = 
            match GeodesicOrder with
            | 3 ->
                [|
                    3.0; -8.0; 16.0; // C1[1]/eps^1: polynomial in eps2 of order 1
                    -1.0; 16.0; // C1[2]/eps^2: polynomial in eps2 of order 0
                    -1.0; 48.0 // C1[3]/eps^3: polynomial in eps2 of order 0
                |]
            | 4 ->
                [|
                    3.0; -8.0; 16.0; // C1[1]/eps^1: polynomial in eps2 of order 1
                    1.0; -2.0; 32.0; // C1[2]/eps^2: polynomial in eps2 of order 1
                    -1.0; 48.0; // C1[3]/eps^3: polynomial in eps2 of order 0
                    -5.0; 512.0 // C1[4]/eps^4: polynomial in eps2 of order 0
                |]
            | 5 -> 
                [|
                    -1.0; 6.0; -16.0; 32.0; // C1[1]/eps^1: polynomial in eps2 of order 2
                    1.0; -2.0; 32.0; // C1[2]/eps^2: polynomial in eps2 of order 1
                    9.0; -16.0; 768.0; // C1[3]/eps^3: polynomial in eps2 of order 1
                    -5.0; 512.0; // C1[4]/eps^4: polynomial in eps2 of order 0
                    -7.0; 1280.0 // C1[5]/eps^5: polynomial in eps2 of order 0
                |]
            | 6 -> 
                [|
                    -1.0; 6.0; -16.0; 32.0; // C1[1]/eps^1: polynomial in eps2 of order 2
                    -9.0; 64.0; -128.0; 2048.0; // C1[2]/eps^2: polynomial in eps2 of order 2
                    9.0; -16.0; 768.0; // C1[3]/eps^3: polynomial in eps2 of order 1
                    3.0; -5.0; 512.0; // C1[4]/eps^4: polynomial in eps2 of order 1
                    -7.0; 1280.0; // C1[5]/eps^5: polynomial in eps2 of order 0
                    -7.0; 2048.0 // C1[6]/eps^6: polynomial in eps2 of order 0
                |]
            | 7 -> 
                [|
                    19.0; -64.0; 384.0; -1024.0; 2048.0; // C1[1]/eps^1: polynomial in eps2 of order 3
                    -9.0; 64.0; -128.0; 2048.0; // C1[2]/eps^2: polynomial in eps2 of order 2
                    -9.0; 72.0; -128.0; 6144.0; // C1[3]/eps^3: polynomial in eps2 of order 2
                    3.0; -5.0; 512.0; // C1[4]/eps^4: polynomial in eps2 of order 1
                    35.0; -56.0; 10240.0; // C1[5]/eps^5: polynomial in eps2 of order 1
                    -7.0; 2048.0; // C1[6]/eps^6: polynomial in eps2 of order 0
                    -33.0; 14336.0 // C1[7]/eps^7: polynomial in eps2 of order 0
                |]
            | 8 -> 
                [|
                    19.0; -64.0; 384.0; -1024.0; 2048.0; // C1[1]/eps^1: polynomial in eps2 of order 3
                    7.0; -18.0; 128.0; -256.0; 4096.0; // C1[2]/eps^2: polynomial in eps2 of order 3
                    -9.0; 72.0; -128.0; 6144.0; // C1[3]/eps^3: polynomial in eps2 of order 2
                    -11.0; 96.0; -160.0; 16384.0; // C1[4]/eps^4: polynomial in eps2 of order 2
                    35.0; -56.0; 10240.0; // C1[5]/eps^5: polynomial in eps2 of order 1
                    9.0; -14.0; 4096.0; // C1[6]/eps^6: polynomial in eps2 of order 1
                    -33.0; 14336.0; // C1[7]/eps^7: polynomial in eps2 of order 0
                    -429.0; 262144.0 // C1[8]/eps^8: polynomial in eps2 of order 0
                |]
            | _ -> raise <| new ArgumentException()
        Fourier nC1 coeff eps

    let C2Fourier eps =
        let coeff = 
            match GeodesicOrder with
            | 3 ->
                [|
                    1.0; 8.0; 16.0; // C2[1]/eps^1: polynomial in eps2 of order 1
                    3.0; 16.0; // C2[2]/eps^2: polynomial in eps2 of order 0
                    5.0; 48.0 // C2[3]/eps^3: polynomial in eps2 of order 0
                |]
            | 4 ->
                [|
                    1.0; 8.0; 16.0; // C2[1]/eps^1.0: polynomial in eps2 of order 1
                    1.0; 6.0; 32.0; // C2[2]/eps^2.0: polynomial in eps2 of order 1
                    5.0; 48.0; // C2[3]/eps^3.0: polynomial in eps2 of order 0
                    35.0; 512.0 // C2[4]/eps^4.0: polynomial in eps2 of order 0
                |]
            | 5 -> 
                [|
                    1.0; 2.0; 16.0; 32.0; // C2[1]/eps^1: polynomial in eps2 of order 2
                    1.0; 6.0; 32.0; // C2[2]/eps^2: polynomial in eps2 of order 1
                    15.0; 80.0; 768.0; // C2[3]/eps^3: polynomial in eps2 of order 1
                    35.0; 512.0; // C2[4]/eps^4: polynomial in eps2 of order 0
                    63.0; 1280.0 // C2[5]/eps^5: polynomial in eps2 of order 0
                |]
            | 6 -> 
                [|
                    1.0; 2.0; 16.0; 32.0; // C2[1]/eps^1: polynomial in eps2 of order 2
                    35.0; 64.0; 384.0; 2048.0; // C2[2]/eps^2: polynomial in eps2 of order 2
                    15.0; 80.0; 768.0; // C2[3]/eps^3: polynomial in eps2 of order 1
                    7.0; 35.0; 512.0; // C2[4]/eps^4: polynomial in eps2 of order 1
                    63.0; 1280.0; // C2[5]/eps^5: polynomial in eps2 of order 0
                    77.0; 2048.0 // C2[6]/eps^6: polynomial in eps2 of order 0
                |]
            | 7 -> 
                [|
                    41.0; 64.0; 128.0; 1024.0; 2048.0; // C2[1]/eps^1: polynomial in eps2 of order 3
                    35.0; 64.0; 384.0; 2048.0; // C2[2]/eps^2: polynomial in eps2 of order 2
                    69.0; 120.0; 640.0; 6144.0; // C2[3]/eps^3: polynomial in eps2 of order 2
                    7.0; 35.0; 512.0; // C2[4]/eps^4: polynomial in eps2 of order 1
                    105.0; 504.0; 10240.0; // C2[5]/eps^5: polynomial in eps2 of order 1
                    77.0; 2048.0; // C2[6]/eps^6: polynomial in eps2 of order 0
                    429.0; 14336.0 // C2[7]/eps^7: polynomial in eps2 of order 0
                |]
            | 8 -> 
                [|
                    41.0; 64.0; 128.0; 1024.0; 2048.0; // C2[1]/eps^1: polynomial in eps2 of order 3
                    47.0; 70.0; 128.0; 768.0; 4096.0; // C2[2]/eps^2: polynomial in eps2 of order 3
                    69.0; 120.0; 640.0; 6144.0; // C2[3]/eps^3: polynomial in eps2 of order 2
                    133.0; 224.0; 1120.0; 16384.0; // C2[4]/eps^4: polynomial in eps2 of order 2
                    105.0; 504.0; 10240.0; // C2[5]/eps^5: polynomial in eps2 of order 1
                    33.0; 154.0; 4096.0; // C2[6]/eps^6: polynomial in eps2 of order 1
                    429.0; 14336.0; // C2[7]/eps^7: polynomial in eps2 of order 0
                    6435.0; 262144.0 // C2[8]/eps^8: polynomial in eps2 of order 0
                |]
            | _ -> raise <| new ArgumentException()
        Fourier nC2 coeff eps

    let C4Coeff =
        match GeodesicOrder with
        | 3 ->
            [|
                -2.0; 105.0; // C4[0]: coeff of eps^2 -- polynomial in n of order 0
                16.0; -7.0; 35.0; // C4[0]: coeff of eps^1 -- polynomial in n of order 1
                8.0; -28.0; 70.0; 105.0; // C4[0]: coeff of eps^0 -- polynomial in n of order 2
                -2.0; 105.0; // C4[1]: coeff of eps^2 -- polynomial in n of order 0
                -16.0; 7.0; 315.0; // C4[1]: coeff of eps^1 -- polynomial in n of order 1
                4.0; 525.0 // C4[2]: coeff of eps^2 -- polynomial in n of order 0
            |]
        | 4 ->
            [|
                11.0; 315.0; // C4[0]: coeff of eps^3 -- polynomial in n of order 0
                -32.0; -6.0; 315.0; // C4[0]: coeff of eps^2 -- polynomial in n of order 1
                -32.0; 48.0; -21.0; 105.0; // C4[0]: coeff of eps^1 -- polynomial in n of order 2
                4.0; 24.0; -84.0; 210.0; 315.0; // C4[0]: coeff of eps^0 -- polynomial in n of order 3
                -1.0; 105.0; // C4[1]: coeff of eps^3 -- polynomial in n of order 0
                64.0; -18.0; 945.0; // C4[1]: coeff of eps^2 -- polynomial in n of order 1
                32.0; -48.0; 21.0; 945.0; // C4[1]: coeff of eps^1 -- polynomial in n of order 2
                -8.0; 1575.0; // C4[2]: coeff of eps^3 -- polynomial in n of order 0
                -32.0; 12.0; 1575.0; // C4[2]: coeff of eps^2 -- polynomial in n of order 1
                8.0; 2205.0 // C4[3]: coeff of eps^3 -- polynomial in n of order 0
            |]
        | 5 ->
            [|
                4.0; 1155.0; // C4[0]: coeff of eps^4 -- polynomial in n of order 0
                -368.0; 121.0; 3465.0; // C4[0]: coeff of eps^3 -- polynomial in n of order 1
                1088.0; -352.0; -66.0; 3465.0; // C4[0]: coeff of eps^2 -- polynomial in n of order 2
                48.0; -352.0; 528.0; -231.0; 1155.0; // C4[0]: coeff of eps^1 -- polynomial in n of order 3
                16.0; 44.0; 264.0; -924.0; 2310.0; 3465.0; // C4[0]: coeff of eps^0 -- polynomial in n of order 4
                4.0; 1155.0; // C4[1]: coeff of eps^4 -- polynomial in n of order 0
                80.0; -99.0; 10395.0; // C4[1]: coeff of eps^3 -- polynomial in n of order 1
                -896.0; 704.0; -198.0; 10395.0; // C4[1]: coeff of eps^2 -- polynomial in n of order 2
                -48.0; 352.0; -528.0; 231.0; 10395.0; // C4[1]: coeff of eps^1 -- polynomial in n of order 3
                -8.0; 1925.0; // C4[2]: coeff of eps^4 -- polynomial in n of order 0
                384.0; -88.0; 17325.0; // C4[2]: coeff of eps^3 -- polynomial in n of order 1
                320.0; -352.0; 132.0; 17325.0; // C4[2]: coeff of eps^2 -- polynomial in n of order 2
                -16.0; 8085.0; // C4[3]: coeff of eps^4 -- polynomial in n of order 0
                -256.0; 88.0; 24255.0; // C4[3]: coeff of eps^3 -- polynomial in n of order 1
                64.0; 31185.0 // C4[4]: coeff of eps^4 -- polynomial in n of order 0
            |]
        | 6 ->
            [|
                97.0; 15015.0; // C4[0]: coeff of eps^5 -- polynomial in n of order 0
                1088.0; 156.0; 45045.0; // C4[0]: coeff of eps^4 -- polynomial in n of order 1
                -224.0; -4784.0; 1573.0; 45045.0; // C4[0]: coeff of eps^3 -- polynomial in n of order 2
                -10656.0; 14144.0; -4576.0; -858.0; 45045.0; // C4[0]: coeff of eps^2 -- polynomial in n of order 3
                64.0; 624.0; -4576.0; 6864.0; -3003.0; 15015.0; // C4[0]: coeff of eps^1 -- polynomial in n of order 4
                100.0; 208.0; 572.0; 3432.0; -12012.0; 30030.0; 45045.0; // C4[0]: coeff of eps^0 -- polynomial in n of order 5
                1.0; 9009.0; // C4[1]: coeff of eps^5 -- polynomial in n of order 0
                -2944.0; 468.0; 135135.0; // C4[1]: coeff of eps^4 -- polynomial in n of order 1
                5792.0; 1040.0; -1287.0; 135135.0; // C4[1]: coeff of eps^3 -- polynomial in n of order 2
                5952.0; -11648.0; 9152.0; -2574.0; 135135.0; // C4[1]: coeff of eps^2 -- polynomial in n of order 3
                -64.0; -624.0; 4576.0; -6864.0; 3003.0; 135135.0; // C4[1]: coeff of eps^1 -- polynomial in n of order 4
                8.0; 10725.0; // C4[2]: coeff of eps^5 -- polynomial in n of order 0
                1856.0; -936.0; 225225.0; // C4[2]: coeff of eps^4 -- polynomial in n of order 1
                -8448.0; 4992.0; -1144.0; 225225.0; // C4[2]: coeff of eps^3 -- polynomial in n of order 2
                -1440.0; 4160.0; -4576.0; 1716.0; 225225.0; // C4[2]: coeff of eps^2 -- polynomial in n of order 3
                -136.0; 63063.0; // C4[3]: coeff of eps^5 -- polynomial in n of order 0
                1024.0; -208.0; 105105.0; // C4[3]: coeff of eps^4 -- polynomial in n of order 1
                3584.0; -3328.0; 1144.0; 315315.0; // C4[3]: coeff of eps^3 -- polynomial in n of order 2
                -128.0; 135135.0; // C4[4]: coeff of eps^5 -- polynomial in n of order 0
                -2560.0; 832.0; 405405.0; // C4[4]: coeff of eps^4 -- polynomial in n of order 1
                128.0; 99099.0 // C4[5]: coeff of eps^5 -- polynomial in n of order 0
            |]
        | 7 ->
            [|
                10.0; 9009.0; // C4[0]: coeff of eps^6 -- polynomial in n of order 0
                -464.0; 291.0; 45045.0; // C4[0]: coeff of eps^5 -- polynomial in n of order 1
                -4480.0; 1088.0; 156.0; 45045.0; // C4[0]: coeff of eps^4 -- polynomial in n of order 2
                10736.0; -224.0; -4784.0; 1573.0; 45045.0; // C4[0]: coeff of eps^3 -- polynomial in n of order 3
                1664.0; -10656.0; 14144.0; -4576.0; -858.0; 45045.0; // C4[0]: coeff of eps^2 -- polynomial in n of order 4
                16.0; 64.0; 624.0; -4576.0; 6864.0; -3003.0; 15015.0; // C4[0]: coeff of eps^1 -- polynomial in n of order 5
                56.0; 100.0; 208.0; 572.0; 3432.0; -12012.0; 30030.0; 45045.0; // C4[0]: coeff of eps^0 -- polynomial in n of order 6
                10.0; 9009.0; // C4[1]: coeff of eps^6 -- polynomial in n of order 0
                112.0; 15.0; 135135.0; // C4[1]: coeff of eps^5 -- polynomial in n of order 1
                3840.0; -2944.0; 468.0; 135135.0; // C4[1]: coeff of eps^4 -- polynomial in n of order 2
                -10704.0; 5792.0; 1040.0; -1287.0; 135135.0; // C4[1]: coeff of eps^3 -- polynomial in n of order 3
                -768.0; 5952.0; -11648.0; 9152.0; -2574.0; 135135.0; // C4[1]: coeff of eps^2 -- polynomial in n of order 4
                -16.0; -64.0; -624.0; 4576.0; -6864.0; 3003.0; 135135.0; // C4[1]: coeff of eps^1 -- polynomial in n of order 5
                -4.0; 25025.0; // C4[2]: coeff of eps^6 -- polynomial in n of order 0
                -1664.0; 168.0; 225225.0; // C4[2]: coeff of eps^5 -- polynomial in n of order 1
                1664.0; 1856.0; -936.0; 225225.0; // C4[2]: coeff of eps^4 -- polynomial in n of order 2
                6784.0; -8448.0; 4992.0; -1144.0; 225225.0; // C4[2]: coeff of eps^3 -- polynomial in n of order 3
                128.0; -1440.0; 4160.0; -4576.0; 1716.0; 225225.0; // C4[2]: coeff of eps^2 -- polynomial in n of order 4
                64.0; 315315.0; // C4[3]: coeff of eps^6 -- polynomial in n of order 0
                1792.0; -680.0; 315315.0; // C4[3]: coeff of eps^5 -- polynomial in n of order 1
                -2048.0; 1024.0; -208.0; 105105.0; // C4[3]: coeff of eps^4 -- polynomial in n of order 2
                -1792.0; 3584.0; -3328.0; 1144.0; 315315.0; // C4[3]: coeff of eps^3 -- polynomial in n of order 3
                -512.0; 405405.0; // C4[4]: coeff of eps^6 -- polynomial in n of order 0
                2048.0; -384.0; 405405.0; // C4[4]: coeff of eps^5 -- polynomial in n of order 1
                3072.0; -2560.0; 832.0; 405405.0; // C4[4]: coeff of eps^4 -- polynomial in n of order 2
                -256.0; 495495.0; // C4[5]: coeff of eps^6 -- polynomial in n of order 0
                -2048.0; 640.0; 495495.0; // C4[5]: coeff of eps^5 -- polynomial in n of order 1
                512.0; 585585.0 // C4[6]: coeff of eps^6 -- polynomial in n of order 0
            |]
        | 8 ->
            [|
                193.0; 85085.0; // C4[0]: coeff of eps^7 -- polynomial in n of order 0
                4192.0; 850.0; 765765.0; // C4[0]: coeff of eps^6 -- polynomial in n of order 1
                20960.0; -7888.0; 4947.0; 765765.0; // C4[0]: coeff of eps^5 -- polynomial in n of order 2
                12480.0; -76160.0; 18496.0; 2652.0; 765765.0; // C4[0]: coeff of eps^4 -- polynomial in n of order 3
                -154048.0; 182512.0; -3808.0; -81328.0; 26741.0; 765765.0; // C4[0]: coeff of eps^3 -- polynomial in n of order 4
                3232.0; 28288.0; -181152.0; 240448.0; -77792.0; -14586.0; 765765.0; // C4[0]: coeff of eps^2 -- polynomial in n of order 5
                96.0; 272.0; 1088.0; 10608.0; -77792.0; 116688.0; -51051.0; 255255.0; // C4[0]: coeff of eps^1 -- polynomial in n of order 6
                588.0; 952.0; 1700.0; 3536.0; 9724.0; 58344.0; -204204.0; 510510.0; 765765.0; // C4[0]: coeff of eps^0 -- polynomial in n of order 7
                349.0; 2297295.0; // C4[1]: coeff of eps^7 -- polynomial in n of order 0
                -1472.0; 510.0; 459459.0; // C4[1]: coeff of eps^6 -- polynomial in n of order 1
                -39840.0; 1904.0; 255.0; 2297295.0; // C4[1]: coeff of eps^5 -- polynomial in n of order 2
                52608.0; 65280.0; -50048.0; 7956.0; 2297295.0; // C4[1]: coeff of eps^4 -- polynomial in n of order 3
                103744.0; -181968.0; 98464.0; 17680.0; -21879.0; 2297295.0; // C4[1]: coeff of eps^3 -- polynomial in n of order 4
                -1344.0; -13056.0; 101184.0; -198016.0; 155584.0; -43758.0; 2297295.0; // C4[1]: coeff of eps^2 -- polynomial in n of order 5
                -96.0; -272.0; -1088.0; -10608.0; 77792.0; -116688.0; 51051.0; 2297295.0; // C4[1]: coeff of eps^1 -- polynomial in n of order 6
                464.0; 1276275.0; // C4[2]: coeff of eps^7 -- polynomial in n of order 0
                -928.0; -612.0; 3828825.0; // C4[2]: coeff of eps^6 -- polynomial in n of order 1
                64256.0; -28288.0; 2856.0; 3828825.0; // C4[2]: coeff of eps^5 -- polynomial in n of order 2
                -126528.0; 28288.0; 31552.0; -15912.0; 3828825.0; // C4[2]: coeff of eps^4 -- polynomial in n of order 3
                -41472.0; 115328.0; -143616.0; 84864.0; -19448.0; 3828825.0; // C4[2]: coeff of eps^3 -- polynomial in n of order 4
                160.0; 2176.0; -24480.0; 70720.0; -77792.0; 29172.0; 3828825.0; // C4[2]: coeff of eps^2 -- polynomial in n of order 5
                -16.0; 97461.0; // C4[3]: coeff of eps^7 -- polynomial in n of order 0
                -16384.0; 1088.0; 5360355.0; // C4[3]: coeff of eps^6 -- polynomial in n of order 1
                -2560.0; 30464.0; -11560.0; 5360355.0; // C4[3]: coeff of eps^5 -- polynomial in n of order 2
                35840.0; -34816.0; 17408.0; -3536.0; 1786785.0; // C4[3]: coeff of eps^4 -- polynomial in n of order 3
                7168.0; -30464.0; 60928.0; -56576.0; 19448.0; 5360355.0; // C4[3]: coeff of eps^3 -- polynomial in n of order 4
                128.0; 2297295.0; // C4[4]: coeff of eps^7 -- polynomial in n of order 0
                26624.0; -8704.0; 6891885.0; // C4[4]: coeff of eps^6 -- polynomial in n of order 1
                -77824.0; 34816.0; -6528.0; 6891885.0; // C4[4]: coeff of eps^5 -- polynomial in n of order 2
                -32256.0; 52224.0; -43520.0; 14144.0; 6891885.0; // C4[4]: coeff of eps^4 -- polynomial in n of order 3
                -6784.0; 8423415.0; // C4[5]: coeff of eps^7 -- polynomial in n of order 0
                24576.0; -4352.0; 8423415.0; // C4[5]: coeff of eps^6 -- polynomial in n of order 1
                45056.0; -34816.0; 10880.0; 8423415.0; // C4[5]: coeff of eps^5 -- polynomial in n of order 2
                -1024.0; 3318315.0; // C4[6]: coeff of eps^7 -- polynomial in n of order 0
                -28672.0; 8704.0; 9954945.0; // C4[6]: coeff of eps^6 -- polynomial in n of order 1
                1024.0; 1640925.0 // C4[7]: coeff of eps^7 -- polynomial in n of order 0
            |]
        | _ -> raise <| new ArgumentException()

    let C3Coeff =
        match GeodesicOrder with
        | 3 ->
            [|
                1.0; 8.0; // C3[1]: coeff of eps^2 -- polynomial in n of order 0
                -1.0; 1.0; 4.0; // C3[1]: coeff of eps^1 -- polynomial in n of order 1
                1.0; 16.0 // C3[2]: coeff of eps^2 -- polynomial in n of order 0
            |]
        | 4 ->
            [|
                3.0; 64.0; // C3[1]: coeff of eps^3 -- polynomial in n of order 0
                // This is a case where a leading 0 term has been inserted to maintain the
                // pattern in the orders of the polynomials.
                0.0; 1.0; 8.0; // C3[1]: coeff of eps^2 -- polynomial in n of order 1
                -1.0; 1.0; 4.0; // C3[1]: coeff of eps^1 -- polynomial in n of order 1
                3.0; 64.0; // C3[2]: coeff of eps^3 -- polynomial in n of order 0
                -3.0; 2.0; 32.0; // C3[2]: coeff of eps^2 -- polynomial in n of order 1
                5.0; 192.0 // C3[3]: coeff of eps^3 -- polynomial in n of order 0
            |]
        | 5 ->
            [|
                5.0; 128.0; // C3[1]: coeff of eps^4 -- polynomial in n of order 0
                3.0; 3.0; 64.0; // C3[1]: coeff of eps^3 -- polynomial in n of order 1
                -1.0; 0.0; 1.0; 8.0; // C3[1]: coeff of eps^2 -- polynomial in n of order 2
                -1.0; 1.0; 4.0; // C3[1]: coeff of eps^1 -- polynomial in n of order 1
                3.0; 128.0; // C3[2]: coeff of eps^4 -- polynomial in n of order 0
                -2.0; 3.0; 64.0; // C3[2]: coeff of eps^3 -- polynomial in n of order 1
                1.0; -3.0; 2.0; 32.0; // C3[2]: coeff of eps^2 -- polynomial in n of order 2
                3.0; 128.0; // C3[3]: coeff of eps^4 -- polynomial in n of order 0
                -9.0; 5.0; 192.0; // C3[3]: coeff of eps^3 -- polynomial in n of order 1
                7.0; 512.0 // C3[4]: coeff of eps^4 -- polynomial in n of order 0
            |]
        | 6->
            [|
                3.0; 128.0; // C3[1]: coeff of eps^5 -- polynomial in n of order 0
                2.0; 5.0; 128.0; // C3[1]: coeff of eps^4 -- polynomial in n of order 1
                -1.0; 3.0; 3.0; 64.0; // C3[1]: coeff of eps^3 -- polynomial in n of order 2
                -1.0; 0.0; 1.0; 8.0; // C3[1]: coeff of eps^2 -- polynomial in n of order 2
                -1.0; 1.0; 4.0; // C3[1]: coeff of eps^1 -- polynomial in n of order 1
                5.0; 256.0; // C3[2]: coeff of eps^5 -- polynomial in n of order 0
                1.0; 3.0; 128.0; // C3[2]: coeff of eps^4 -- polynomial in n of order 1
                -3.0; -2.0; 3.0; 64.0; // C3[2]: coeff of eps^3 -- polynomial in n of order 2
                1.0; -3.0; 2.0; 32.0; // C3[2]: coeff of eps^2 -- polynomial in n of order 2
                7.0; 512.0; // C3[3]: coeff of eps^5 -- polynomial in n of order 0
                -10.0; 9.0; 384.0; // C3[3]: coeff of eps^4 -- polynomial in n of order 1
                5.0; -9.0; 5.0; 192.0; // C3[3]: coeff of eps^3 -- polynomial in n of order 2
                7.0; 512.0; // C3[4]: coeff of eps^5 -- polynomial in n of order 0
                -14.0; 7.0; 512.0; // C3[4]: coeff of eps^4 -- polynomial in n of order 1
                21.0; 2560.0 // C3[5]: coeff of eps^5 -- polynomial in n of order 0
            |]
        | 7->
            [|
                21.0; 1024.0; // C3[1]: coeff of eps^6 -- polynomial in n of order 0
                11.0; 12.0; 512.0; // C3[1]: coeff of eps^5 -- polynomial in n of order 1
                2.0; 2.0; 5.0; 128.0; // C3[1]: coeff of eps^4 -- polynomial in n of order 2
                -5.0; -1.0; 3.0; 3.0; 64.0; // C3[1]: coeff of eps^3 -- polynomial in n of order 3
                -1.0; 0.0; 1.0; 8.0; // C3[1]: coeff of eps^2 -- polynomial in n of order 2
                -1.0; 1.0; 4.0; // C3[1]: coeff of eps^1 -- polynomial in n of order 1
                27.0; 2048.0; // C3[2]: coeff of eps^6 -- polynomial in n of order 0
                1.0; 5.0; 256.0; // C3[2]: coeff of eps^5 -- polynomial in n of order 1
                -9.0; 2.0; 6.0; 256.0; // C3[2]: coeff of eps^4 -- polynomial in n of order 2
                2.0; -3.0; -2.0; 3.0; 64.0; // C3[2]: coeff of eps^3 -- polynomial in n of order 3
                1.0; -3.0; 2.0; 32.0; // C3[2]: coeff of eps^2 -- polynomial in n of order 2
                3.0; 256.0; // C3[3]: coeff of eps^6 -- polynomial in n of order 0
                -4.0; 21.0; 1536.0; // C3[3]: coeff of eps^5 -- polynomial in n of order 1
                -6.0; -10.0; 9.0; 384.0; // C3[3]: coeff of eps^4 -- polynomial in n of order 2
                -1.0; 5.0; -9.0; 5.0; 192.0; // C3[3]: coeff of eps^3 -- polynomial in n of order 3
                9.0; 1024.0; // C3[4]: coeff of eps^6 -- polynomial in n of order 0
                -10.0; 7.0; 512.0; // C3[4]: coeff of eps^5 -- polynomial in n of order 1
                10.0; -14.0; 7.0; 512.0; // C3[4]: coeff of eps^4 -- polynomial in n of order 2
                9.0; 1024.0; // C3[5]: coeff of eps^6 -- polynomial in n of order 0
                -45.0; 21.0; 2560.0; // C3[5]: coeff of eps^5 -- polynomial in n of order 1
                11.0; 2048.0; // C3[6]: coeff of eps^6 -- polynomial in n of order 0
            |]
        | 8->
            [|
                243.0; 16384.0; // C3[1]: coeff of eps^7 -- polynomial in n of order 0
                10.0; 21.0; 1024.0; // C3[1]: coeff of eps^6 -- polynomial in n of order 1
                3.0; 11.0; 12.0; 512.0; // C3[1]: coeff of eps^5 -- polynomial in n of order 2
                -2.0; 2.0; 2.0; 5.0; 128.0; // C3[1]: coeff of eps^4 -- polynomial in n of order 3
                -5.0; -1.0; 3.0; 3.0; 64.0; // C3[1]: coeff of eps^3 -- polynomial in n of order 3
                -1.0; 0.0; 1.0; 8.0; // C3[1]: coeff of eps^2 -- polynomial in n of order 2
                -1.0; 1.0; 4.0; // C3[1]: coeff of eps^1 -- polynomial in n of order 1
                187.0; 16384.0; // C3[2]: coeff of eps^7 -- polynomial in n of order 0
                69.0; 108.0; 8192.0; // C3[2]: coeff of eps^6 -- polynomial in n of order 1
                -2.0; 1.0; 5.0; 256.0; // C3[2]: coeff of eps^5 -- polynomial in n of order 2
                -6.0; -9.0; 2.0; 6.0; 256.0; // C3[2]: coeff of eps^4 -- polynomial in n of order 3
                2.0; -3.0; -2.0; 3.0; 64.0; // C3[2]: coeff of eps^3 -- polynomial in n of order 3
                1.0; -3.0; 2.0; 32.0; // C3[2]: coeff of eps^2 -- polynomial in n of order 2
                139.0; 16384.0; // C3[3]: coeff of eps^7 -- polynomial in n of order 0
                -1.0; 12.0; 1024.0; // C3[3]: coeff of eps^6 -- polynomial in n of order 1
                -77.0; -8.0; 42.0; 3072.0; // C3[3]: coeff of eps^5 -- polynomial in n of order 2
                10.0; -6.0; -10.0; 9.0; 384.0; // C3[3]: coeff of eps^4 -- polynomial in n of order 3
                -1.0; 5.0; -9.0; 5.0; 192.0; // C3[3]: coeff of eps^3 -- polynomial in n of order 3
                127.0; 16384.0; // C3[4]: coeff of eps^7 -- polynomial in n of order 0
                -43.0; 72.0; 8192.0; // C3[4]: coeff of eps^6 -- polynomial in n of order 1
                -7.0; -40.0; 28.0; 2048.0; // C3[4]: coeff of eps^5 -- polynomial in n of order 2
                -7.0; 20.0; -28.0; 14.0; 1024.0; // C3[4]: coeff of eps^4 -- polynomial in n of order 3
                99.0; 16384.0; // C3[5]: coeff of eps^7 -- polynomial in n of order 0
                -15.0; 9.0; 1024.0; // C3[5]: coeff of eps^6 -- polynomial in n of order 1
                75.0; -90.0; 42.0; 5120.0; // C3[5]: coeff of eps^5 -- polynomial in n of order 2
                99.0; 16384.0; // C3[6]: coeff of eps^7 -- polynomial in n of order 0
                -99.0; 44.0; 8192.0; // C3[6]: coeff of eps^6 -- polynomial in n of order 1
                429.0; 114688.0; // C3[7]: coeff of eps^7 -- polynomial in n of order 0
            |]
        | _ -> raise <| new ArgumentException()

    let A3Coeff =
        match GeodesicOrder with
        | 3 -> 
            [| 
                -1.0; 4.0; // A3 coeff of eps^2
                1.0; -1.0; 2.0; // A3 coeff of eps^1
                1.0; 1.0 // A3 coeff of eps^0
            |]
        | 4 -> 
            [|
                -1.0; 16.0; // A3 coeff of eps^3
                -1.0; -2.0; 8.0; // A3 coeff of eps^2
                1.0; -1.0; 2.0; // A3 coeff of eps^1
                1.0; 1.0 // A3 coeff of eps^0
            |]
        | 5 -> 
            [|
                -3.0; 64.0; // A3 coeff of eps^4
                -3.0; -1.0; 16.0; // A3 coeff of eps^3
                3.0; -1.0; -2.0; 8.0; // A3 coeff of eps^2
                1.0; -1.0; 2.0; // A3 coeff of eps^1
                1.0; 1.0 // A3 coeff of eps^0
            |]
        | 6 -> 
            [|
                -3.0; 128.0; // A3 coeff of eps^5
                -2.0; -3.0; 64.0; // A3 coeff of eps^4
                -1.0; -3.0; -1.0; 16.0; // A3 coeff of eps^3
                3.0; -1.0; -2.0; 8.0; // A3 coeff of eps^2
                1.0; -1.0; 2.0; // A3 coeff of eps^1
                1.0; 1.0 // A3 coeff of eps^0
            |]
        | 7 -> 
            [|
                -5.0; 256.0; // A3 coeff of eps^6
                -5.0; -3.0; 128.0; // A3 coeff of eps^5
                -10.0; -2.0; -3.0; 64.0; // A3 coeff of eps^4
                5.0; -1.0; -3.0; -1.0; 16.0; // A3 coeff of eps^3
                3.0; -1.0; -2.0; 8.0; // A3 coeff of eps^2
                1.0; -1.0; 2.0; // A3 coeff of eps^1
                1.0; 1.0 // A3 coeff of eps^0
            |]
        | 8 -> 
            [|
                -25.0; 2048.0; // A3 coeff of eps^7
                -15.0; -20.0; 1024.0; // A3 coeff of eps^6
                -5.0; -10.0; -6.0; 256.0; // A3 coeff of eps^5
                -5.0; -20.0; -4.0; -6.0; 128.0; // A3 coeff of eps^4
                5.0; -1.0; -3.0; -1.0; 16.0; // A3 coeff of eps^3
                3.0; -1.0; -2.0; 8.0; // A3 coeff of eps^2
                1.0; -1.0; 2.0; // A3 coeff of eps^1
                1.0; 1.0 // A3 coeff of eps^0
            |]
        | _ -> raise <| new ArgumentException()

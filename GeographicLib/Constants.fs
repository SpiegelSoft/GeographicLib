namespace GeographicLib

type LowToHighRatio =
    struct
        val Ratio : float
        new(ratio : float) = { Ratio = if ratio <= 1.0 then ratio else 1.0 / ratio }
    end

module Constants = 
    [<Literal>]
    let WGS84_a = 6378137.0<m>
    let WGS84_f = LowToHighRatio(1.0 / ( 298257223563.0 / 1000000000.0 ))
    // The precision of floating point numbers used in %GeographicLib.  1 means
    // float (single precision); 2 (the default) means double; 3 means long double;
    // 4 is reserved for quadruple precision.  Nearly all the testing has been
    // carried out with doubles and that's the recommended configuration.  In order
    // for long double to be used, GEOGRAPHICLIB_HAVE_LONG_DOUBLE needs to be
    // defined.  Note that with Microsoft Visual Studio, long double is the same as
    // double.
    [<Literal>]
    let Precision = 2


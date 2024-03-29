using System;
using System.Linq;
using System.Collections.Generic;
using Core;

namespace ProjectC
{
    public static partial class MainClass
    {
        //For info: numpy.random.seed(246754954)



        // data for SquareSubMatrix() test.
        // a size 4 list of arguments for the test routine.
        private static readonly List<object> SquareSubMatrix_L = new List<object>();
        private static readonly double[,] Temp001 =
        {
            {7.190000, -4.370000, 7.220000, -3.380000, 7.740000, -9.070000, 9.890000, 9.020000, -2.800000, -7.690000},
            {1.370000, 3.680000, 3.000000, -0.070000, -5.360000, -9.670000, 9.310000, -1.290000, -7.880000, -10.000000},
            {1.870000, -8.880000, -9.260000, 5.510000, -6.730000, 6.600000, -8.850000, -0.120000, -1.370000, 8.460000},
            {-2.930000, 8.730000, 3.130000, 3.290000, -1.450000, -0.380000, -9.820000, 1.820000, -4.780000, 2.190000},
            {0.450000, -5.800000, -7.090000, -3.410000, -9.370000, -9.090000, -7.150000, -3.570000, 9.200000, -4.000000},
            {-6.280000, -5.790000, -8.730000, -6.520000, 8.260000, 0.030000, -2.790000, -9.950000, 9.080000, 9.870000},
            {6.950000, -2.870000, 5.150000, 0.910000, -2.890000, -0.420000, -6.230000, 3.170000, -9.540000, 4.500000},
            {-3.150000, -6.220000, 2.100000, -0.330000, -7.410000, 0.430000, -5.450000, 0.350000, -7.820000, -8.460000},
            {0.890000, -4.660000, -3.200000, -8.040000, -2.930000, 7.010000, 8.380000, -1.690000, -6.880000, 4.670000},
            {8.810000, 2.480000, -6.470000, 4.900000, -8.070000, 9.730000, -6.480000, -3.610000, -3.930000, -0.920000}
        };
        private static readonly int Temp002 = 2;
        private static readonly int Temp003 = 0;
        private static readonly double[,] Temp004 =
        {
            {-4.370000, 7.220000, -3.380000, 7.740000, -9.070000, 9.890000, 9.020000, -2.800000, -7.690000},
            {3.680000, 3.000000, -0.070000, -5.360000, -9.670000, 9.310000, -1.290000, -7.880000, -10.000000},
            {8.730000, 3.130000, 3.290000, -1.450000, -0.380000, -9.820000, 1.820000, -4.780000, 2.190000},
            {-5.800000, -7.090000, -3.410000, -9.370000, -9.090000, -7.150000, -3.570000, 9.200000, -4.000000},
            {-5.790000, -8.730000, -6.520000, 8.260000, 0.030000, -2.790000, -9.950000, 9.080000, 9.870000},
            {-2.870000, 5.150000, 0.910000, -2.890000, -0.420000, -6.230000, 3.170000, -9.540000, 4.500000},
            {-6.220000, 2.100000, -0.330000, -7.410000, 0.430000, -5.450000, 0.350000, -7.820000, -8.460000},
            {-4.660000, -3.200000, -8.040000, -2.930000, 7.010000, 8.380000, -1.690000, -6.880000, 4.670000},
            {2.480000, -6.470000, 4.900000, -8.070000, 9.730000, -6.480000, -3.610000, -3.930000, -0.920000}
        };
        private static readonly double[,] Temp006 =
        {
            {-0.700000, -2.520000, -9.890000, 2.870000, -0.480000, 3.030000, -7.730000, 2.150000, -3.480000, -8.070000, -8.020000, 8.460000, 9.980000},
            {9.080000, 6.320000, 0.260000, -4.660000, -2.580000, 9.050000, 1.590000, 4.060000, 5.970000, 8.260000, -5.920000, 6.800000, 0.100000},
            {5.450000, 7.140000, -0.190000, -6.430000, 3.430000, -0.760000, 6.560000, 2.720000, 3.690000, 6.280000, 2.820000, 6.750000, 0.690000},
            {9.400000, -6.090000, -7.220000, 9.520000, -0.950000, 6.610000, 8.430000, 2.140000, 5.210000, 2.170000, 5.950000, -9.220000, 3.220000},
            {-0.030000, -1.750000, -4.860000, 1.220000, 1.720000, 1.790000, -0.850000, -4.930000, 1.480000, 6.610000, 8.400000, 4.600000, 4.970000},
            {0.020000, 6.080000, 2.390000, -6.570000, -8.620000, -1.710000, -6.560000, 6.770000, 2.180000, -6.420000, 3.040000, -7.740000, 6.540000},
            {2.570000, -6.920000, 6.820000, -3.410000, -7.960000, -2.200000, 5.330000, 8.750000, 1.910000, 8.710000, -7.980000, -1.980000, 7.710000},
            {0.560000, 8.510000, 4.940000, -8.770000, -7.260000, 3.580000, -3.570000, 6.660000, -3.850000, 4.460000, -1.640000, -8.740000, 4.910000},
            {-2.910000, -6.780000, -5.090000, 0.030000, -2.190000, 4.510000, -1.670000, -1.280000, 8.930000, -9.400000, 7.440000, 8.580000, 5.170000},
            {-1.370000, 1.020000, -2.110000, 9.520000, 9.710000, 6.630000, 1.240000, -2.500000, 9.000000, 3.370000, 6.420000, -1.600000, 9.150000},
            {-8.320000, 7.410000, 4.910000, 5.760000, -6.970000, 1.400000, 4.470000, 8.830000, 6.890000, -6.110000, -7.230000, 1.070000, -7.910000},
            {-9.210000, -1.140000, 2.710000, 2.520000, -4.660000, -4.530000, -2.280000, -4.490000, 6.700000, -2.960000, -8.590000, -4.720000, 9.630000},
            {2.380000, 3.900000, -0.260000, -1.640000, -8.780000, 8.230000, 2.570000, -0.120000, 5.700000, 7.470000, 0.110000, -1.510000, -9.750000}
        };
        private static readonly int Temp007 = 6;
        private static readonly int Temp008 = 3;
        private static readonly double[,] Temp009 =
        {
            {-0.700000, -2.520000, -9.890000, -0.480000, 3.030000, -7.730000, 2.150000, -3.480000, -8.070000, -8.020000, 8.460000, 9.980000},
            {9.080000, 6.320000, 0.260000, -2.580000, 9.050000, 1.590000, 4.060000, 5.970000, 8.260000, -5.920000, 6.800000, 0.100000},
            {5.450000, 7.140000, -0.190000, 3.430000, -0.760000, 6.560000, 2.720000, 3.690000, 6.280000, 2.820000, 6.750000, 0.690000},
            {9.400000, -6.090000, -7.220000, -0.950000, 6.610000, 8.430000, 2.140000, 5.210000, 2.170000, 5.950000, -9.220000, 3.220000},
            {-0.030000, -1.750000, -4.860000, 1.720000, 1.790000, -0.850000, -4.930000, 1.480000, 6.610000, 8.400000, 4.600000, 4.970000},
            {0.020000, 6.080000, 2.390000, -8.620000, -1.710000, -6.560000, 6.770000, 2.180000, -6.420000, 3.040000, -7.740000, 6.540000},
            {0.560000, 8.510000, 4.940000, -7.260000, 3.580000, -3.570000, 6.660000, -3.850000, 4.460000, -1.640000, -8.740000, 4.910000},
            {-2.910000, -6.780000, -5.090000, -2.190000, 4.510000, -1.670000, -1.280000, 8.930000, -9.400000, 7.440000, 8.580000, 5.170000},
            {-1.370000, 1.020000, -2.110000, 9.710000, 6.630000, 1.240000, -2.500000, 9.000000, 3.370000, 6.420000, -1.600000, 9.150000},
            {-8.320000, 7.410000, 4.910000, -6.970000, 1.400000, 4.470000, 8.830000, 6.890000, -6.110000, -7.230000, 1.070000, -7.910000},
            {-9.210000, -1.140000, 2.710000, -4.660000, -4.530000, -2.280000, -4.490000, 6.700000, -2.960000, -8.590000, -4.720000, 9.630000},
            {2.380000, 3.900000, -0.260000, -8.780000, 8.230000, 2.570000, -0.120000, 5.700000, 7.470000, 0.110000, -1.510000, -9.750000}
        };
        private static readonly double[,] Temp011 =
        {
            {8.200000, 0.150000, 3.500000, 6.420000, 8.340000, 8.980000, -6.740000, -9.960000, 6.560000, 2.110000},
            {-3.270000, -6.630000, -3.490000, 1.920000, -6.240000, -9.690000, -5.960000, -8.840000, 0.230000, 4.670000},
            {1.630000, -9.780000, 3.030000, -9.450000, -2.970000, -0.280000, 3.110000, 1.750000, 8.200000, 9.380000},
            {-1.730000, -4.230000, 7.100000, 9.030000, -1.300000, -6.750000, 7.360000, 5.960000, -7.380000, -4.900000},
            {4.680000, -3.930000, 6.230000, 9.560000, -4.670000, -0.300000, 7.160000, 8.890000, 8.520000, 3.740000},
            {3.570000, -1.960000, 6.220000, 0.410000, 3.090000, -9.230000, -9.350000, 9.670000, -1.710000, 1.660000},
            {1.640000, -2.720000, 8.900000, -8.720000, 7.970000, 6.180000, -2.750000, -3.160000, 3.330000, -0.490000},
            {-4.890000, 6.690000, -3.710000, -9.610000, -6.840000, -8.800000, -2.930000, -0.660000, -0.790000, -3.030000},
            {0.200000, 6.830000, 6.270000, -0.320000, -4.670000, -2.040000, -5.870000, 1.850000, -1.670000, 7.390000},
            {-7.270000, 1.520000, 0.810000, 1.930000, -7.080000, -0.270000, -9.460000, 5.250000, 7.970000, 5.100000}
        };
        private static readonly int Temp012 = 9;
        private static readonly int Temp013 = 1;
        private static readonly double[,] Temp014 =
        {
            {8.200000, 3.500000, 6.420000, 8.340000, 8.980000, -6.740000, -9.960000, 6.560000, 2.110000},
            {-3.270000, -3.490000, 1.920000, -6.240000, -9.690000, -5.960000, -8.840000, 0.230000, 4.670000},
            {1.630000, 3.030000, -9.450000, -2.970000, -0.280000, 3.110000, 1.750000, 8.200000, 9.380000},
            {-1.730000, 7.100000, 9.030000, -1.300000, -6.750000, 7.360000, 5.960000, -7.380000, -4.900000},
            {4.680000, 6.230000, 9.560000, -4.670000, -0.300000, 7.160000, 8.890000, 8.520000, 3.740000},
            {3.570000, 6.220000, 0.410000, 3.090000, -9.230000, -9.350000, 9.670000, -1.710000, 1.660000},
            {1.640000, 8.900000, -8.720000, 7.970000, 6.180000, -2.750000, -3.160000, 3.330000, -0.490000},
            {-4.890000, -3.710000, -9.610000, -6.840000, -8.800000, -2.930000, -0.660000, -0.790000, -3.030000},
            {0.200000, 6.270000, -0.320000, -4.670000, -2.040000, -5.870000, 1.850000, -1.670000, 7.390000}
        };
        private static readonly double[,] Temp016 =
        {
            {1.390000, -5.120000, 8.010000, -6.600000, -7.710000, -9.990000, 6.430000, -9.150000, -4.360000, 9.490000},
            {0.480000, -6.810000, -6.680000, 8.370000, -5.100000, -9.270000, 8.220000, -8.870000, -4.160000, -2.140000},
            {-0.980000, -7.210000, 6.650000, 2.290000, 6.890000, -5.790000, 8.810000, -8.740000, -1.900000, -8.730000},
            {-8.110000, -9.530000, 8.500000, -2.810000, 6.000000, -3.960000, 4.110000, -0.640000, 5.740000, -3.830000},
            {-1.950000, -8.330000, -0.180000, -5.640000, 4.850000, -4.110000, -3.050000, 3.440000, 4.220000, -6.330000},
            {-4.040000, 7.040000, 8.740000, -4.530000, -4.570000, 8.440000, 0.960000, -8.600000, 9.680000, -0.500000},
            {6.010000, -7.060000, -2.340000, 3.910000, -1.170000, 9.930000, -8.410000, -5.720000, -9.550000, -2.480000},
            {-6.080000, 6.420000, 9.830000, 5.130000, -7.920000, 3.060000, 9.060000, -3.420000, 6.840000, 5.710000},
            {4.940000, -9.910000, 0.900000, 0.770000, -2.370000, 9.760000, -3.210000, 4.710000, -1.810000, -3.590000},
            {-3.360000, -0.950000, -5.650000, -7.780000, -6.250000, -4.890000, -0.500000, 7.920000, 8.990000, -3.580000}
        };
        private static readonly int Temp017 = 4;
        private static readonly int Temp018 = 3;
        private static readonly double[,] Temp019 =
        {
            {1.390000, -5.120000, 8.010000, -7.710000, -9.990000, 6.430000, -9.150000, -4.360000, 9.490000},
            {0.480000, -6.810000, -6.680000, -5.100000, -9.270000, 8.220000, -8.870000, -4.160000, -2.140000},
            {-0.980000, -7.210000, 6.650000, 6.890000, -5.790000, 8.810000, -8.740000, -1.900000, -8.730000},
            {-8.110000, -9.530000, 8.500000, 6.000000, -3.960000, 4.110000, -0.640000, 5.740000, -3.830000},
            {-4.040000, 7.040000, 8.740000, -4.570000, 8.440000, 0.960000, -8.600000, 9.680000, -0.500000},
            {6.010000, -7.060000, -2.340000, -1.170000, 9.930000, -8.410000, -5.720000, -9.550000, -2.480000},
            {-6.080000, 6.420000, 9.830000, -7.920000, 3.060000, 9.060000, -3.420000, 6.840000, 5.710000},
            {4.940000, -9.910000, 0.900000, -2.370000, 9.760000, -3.210000, 4.710000, -1.810000, -3.590000},
            {-3.360000, -0.950000, -5.650000, -6.250000, -4.890000, -0.500000, 7.920000, 8.990000, -3.580000}
        };


        // data for Determinant() test.
        // a size 4 list of arguments for the test routine.
        private static readonly List<object> Determinant_L = new List<object>();
        private static readonly double[,] Temp021 =
        {
            {-6.140000, 5.760000, -0.720000, 4.630000, 5.430000, -0.090000},
            {-4.310000, 1.510000, 7.840000, 6.530000, 7.640000, 8.070000},
            {3.040000, 7.840000, -0.840000, -7.570000, 5.030000, -6.880000},
            {2.410000, 8.720000, -6.030000, -8.890000, 4.830000, 3.170000},
            {6.030000, 3.080000, -8.450000, -3.110000, -6.940000, 0.390000},
            {2.950000, -0.580000, 5.680000, -2.980000, 1.620000, 8.410000}
        };
        private static readonly double Temp022 = -169205.047336;
        private static readonly double[,] Temp024 =
        {
            {-6.300000, 0.320000, -5.980000, -4.970000, 8.190000},
            {-8.920000, 1.220000, 4.630000, -5.290000, -2.520000},
            {-2.040000, -1.190000, 4.350000, 1.200000, 4.870000},
            {7.670000, 1.600000, -3.810000, 0.530000, 7.110000},
            {-4.850000, -8.550000, 1.560000, 8.840000, 2.650000}
        };
        private static readonly double Temp025 = -3204.174309;
        private static readonly double[,] Temp027 =
        {
            {7.900000, 5.340000, 6.170000, -2.190000},
            {9.930000, -6.030000, -4.030000, 7.890000},
            {-7.480000, 4.150000, 7.860000, 0.120000},
            {-5.250000, 7.580000, -3.810000, 4.370000}
        };
        private static readonly double Temp028 = -12150.712408;
        private static readonly double[,] Temp030 =
        {
            {-2.620000, -8.980000, -8.340000, -1.070000, 1.120000, -3.660000},
            {-3.760000, -3.410000, -5.780000, -2.160000, 6.930000, 3.350000},
            {3.640000, 0.650000, -6.780000, -1.800000, -3.520000, -7.530000},
            {-4.810000, 4.960000, 7.890000, -7.010000, 2.850000, -7.160000},
            {3.450000, -0.970000, 2.880000, 6.450000, -9.620000, -9.700000},
            {8.330000, 3.110000, 7.430000, -0.890000, -7.170000, 9.020000}
        };
        private static readonly double Temp031 = -84410.101197;


        // data for GramSchmidt() test.
        // a size 4 list of arguments for the test routine.
        private static readonly List<object> GramSchmidt_L = new List<object>();
        private static readonly double[,] Temp033 =
        {
            {6.730000, -6.200000, 5.300000, -0.220000, 9.470000},
            {-1.300000, 6.290000, -8.580000, -2.410000, -8.790000},
            {1.140000, 0.790000, -0.950000, -9.000000, 6.860000},
            {7.430000, -9.840000, 8.630000, 8.850000, 3.630000},
            {-3.400000, 0.170000, -7.700000, -8.230000, 0.150000},
            {-0.150000, -7.480000, -2.850000, -9.570000, 6.390000},
            {-5.540000, -1.730000, 4.360000, 0.630000, -8.960000}
        };
        private static readonly double[,] Temp034 =
        {
            {0.557434, -0.084512, -0.033539, -0.315322, -0.173396},
            {-0.107677, 0.434365, -0.417774, 0.236203, -0.742930},
            {0.094424, 0.136751, -0.084467, -0.807862, -0.209731},
            {0.615414, -0.339176, 0.071898, 0.401982, -0.276683},
            {-0.281616, -0.200479, -0.531320, 0.061451, 0.129461},
            {-0.012424, -0.623504, -0.549478, -0.129667, -0.013750},
            {-0.458869, -0.491421, 0.477293, -0.099801, -0.529621}
        };
        private static readonly double[,] Temp035 =
        {
            {12.073173, -9.275557, 9.302774, 6.880975, 13.096914},
            {0.000000, 12.181545, -6.053689, 2.046586, -4.522631},
            {0.000000, 0.000000, 11.845638, 12.342673, -4.831249},
            {0.000000, 0.000000, 0.000000, 11.000718, -9.070195},
            {0.000000, 0.000000, 0.000000, 0.000000, 7.122146}
        };
        private static readonly double[,] Temp037 =
        {
            {-1.670000, 6.290000, -2.030000, 8.640000},
            {5.420000, -3.490000, 9.990000, 4.420000},
            {5.360000, -8.140000, -5.720000, 4.020000},
            {-8.270000, -4.440000, 6.530000, 8.230000}
        };
        private static readonly double[,] Temp038 =
        {
            {-0.146872, 0.515475, -0.051774, 0.842635},
            {0.476673, -0.174187, 0.827412, 0.240480},
            {0.471396, -0.587478, -0.514394, 0.409944},
            {-0.727322, -0.599011, 0.219333, 0.253144}
        };
        private static readonly double[,] Temp039 =
        {
            {11.370479, -3.195265, -2.385687, -3.252924},
            {0.000000, 11.291930, -3.337702, -3.607723},
            {0.000000, 0.000000, 12.745530, 2.947083},
            {0.000000, 0.000000, 0.000000, 12.074636}
        };
        private static readonly double[,] Temp041 =
        {
            {-9.850000, 4.010000, -6.590000, -2.180000, 7.230000},
            {-3.120000, -4.750000, -5.210000, 6.790000, -7.440000},
            {-9.400000, -3.030000, 2.370000, -4.690000, -2.570000},
            {1.940000, -2.130000, -0.110000, -1.810000, -6.550000},
            {-0.640000, -1.710000, 5.630000, -6.570000, 5.490000},
            {-0.350000, -7.200000, 9.070000, -4.910000, 0.130000},
            {-9.960000, -6.460000, 0.760000, 8.500000, 1.750000}
        };
        private static readonly double[,] Temp042 =
        {
            {-0.570008, 0.540044, -0.049661, -0.323170, 0.426340},
            {-0.180551, -0.349898, -0.732003, -0.273862, 0.150284},
            {-0.543967, -0.077984, 0.293641, -0.303967, -0.679740},
            {0.112266, -0.222340, -0.217178, -0.489337, -0.162842},
            {-0.037036, -0.135433, 0.432846, -0.229373, 0.473658},
            {-0.020254, -0.616195, 0.366445, -0.228174, 0.269242},
            {-0.576374, -0.363825, -0.081709, 0.619709, 0.101198}
        };
        private static readonly double[,] Temp043 =
        {
            {17.280457, 3.913513, 2.565210, -2.191719, -3.329819},
            {0.000000, 11.555973, -8.524170, -1.962126, 6.704177},
            {0.000000, 0.000000, 10.559300, -11.183694, 8.035886},
            {0.000000, 0.000000, 0.000000, 9.051139, 3.482933},
            {0.000000, 0.000000, 0.000000, 0.000000, 7.590350}
        };
        private static readonly double[,] Temp045 =
        {
            {-3.860000, -6.540000, -9.770000, 1.220000, -4.830000, -1.990000, 6.060000},
            {4.640000, 2.680000, -8.590000, -5.580000, 0.310000, -8.010000, 3.370000},
            {-5.550000, -8.810000, 4.870000, -4.490000, -5.430000, -8.830000, 7.780000},
            {-2.330000, -0.140000, -4.160000, 1.690000, -2.850000, -2.390000, 7.050000},
            {0.120000, 3.890000, -1.980000, -0.410000, 6.620000, 0.850000, -0.610000},
            {-2.250000, 0.000000, 4.700000, 0.970000, -5.650000, -0.140000, 1.700000},
            {-6.760000, -9.580000, 9.850000, -1.330000, -5.990000, -8.370000, 4.290000},
            {-2.460000, 8.210000, 8.880000, 7.800000, -7.790000, 3.860000, -4.810000},
            {-9.680000, -8.650000, -8.500000, -8.620000, 4.270000, -2.500000, 7.580000},
            {8.950000, 5.730000, 4.570000, 5.780000, 5.100000, -0.110000, -2.190000},
            {-3.380000, -8.450000, -2.140000, 9.320000, -9.550000, 3.620000, 2.720000},
            {8.590000, 1.900000, -7.920000, -7.220000, 3.100000, -7.980000, 3.910000}
        };
        private static readonly double[,] Temp046 =
        {
            {-0.195835, -0.224862, -0.382708, 0.331699, -0.161274, -0.206396, -0.202088},
            {0.235407, -0.068224, -0.328823, -0.181819, -0.420904, -0.313505, -0.337047},
            {-0.281576, -0.284546, 0.277809, -0.271976, -0.218028, -0.199104, 0.268305},
            {-0.118211, 0.113036, -0.238573, 0.171909, -0.143976, -0.460306, 0.540524},
            {0.006088, 0.247851, -0.171794, -0.038827, 0.368106, -0.390554, -0.079104},
            {-0.114152, 0.117987, 0.154274, -0.044157, -0.318165, 0.142177, 0.487753},
            {-0.342964, -0.271401, 0.486992, -0.174187, -0.052930, -0.398003, -0.351231},
            {-0.124807, 0.665380, 0.152157, 0.172031, -0.398862, -0.189483, -0.153993},
            {-0.491109, -0.057520, -0.420827, -0.284811, 0.318931, -0.064877, 0.140664},
            {0.454073, -0.094971, 0.293342, 0.226373, 0.348842, -0.474247, 0.208967},
            {-0.171482, -0.374818, 0.010839, 0.699138, -0.097037, 0.098233, 0.007408},
            {0.435808, -0.326317, -0.185545, -0.262381, -0.311637, -0.042014, 0.163277}
        };
        private static readonly double[,] Temp047 =
        {
            {19.710507, 15.820486, -2.858506, 0.994581, 9.803974, 0.762258, -7.965696},
            {0.000000, 15.306287, 7.667500, 5.948356, 1.542358, 9.671895, -9.848262},
            {0.000000, 0.000000, 22.486220, 7.240283, -6.172914, 0.395868, -5.751840},
            {0.000000, 0.000000, 0.000000, 16.651386, -8.527656, 10.193950, -2.911641},
            {0.000000, 0.000000, 0.000000, 0.000000, 13.002936, 6.522402, -4.009931},
            {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 8.933367, -6.767791},
            {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 4.916125}
        };


        private static void PopulateSquareSubMatrix_L()
        {
            var Temp000 = new List<object>();
            Temp000.Add(new Matrix(Temp001));
            Temp000.Add(Temp002);
            Temp000.Add(Temp003);
            Temp000.Add(new Matrix(Temp004));
            SquareSubMatrix_L.Add(Temp000);
            var Temp005 = new List<object>();
            Temp005.Add(new Matrix(Temp006));
            Temp005.Add(Temp007);
            Temp005.Add(Temp008);
            Temp005.Add(new Matrix(Temp009));
            SquareSubMatrix_L.Add(Temp005);
            var Temp010 = new List<object>();
            Temp010.Add(new Matrix(Temp011));
            Temp010.Add(Temp012);
            Temp010.Add(Temp013);
            Temp010.Add(new Matrix(Temp014));
            SquareSubMatrix_L.Add(Temp010);
            var Temp015 = new List<object>();
            Temp015.Add(new Matrix(Temp016));
            Temp015.Add(Temp017);
            Temp015.Add(Temp018);
            Temp015.Add(new Matrix(Temp019));
            SquareSubMatrix_L.Add(Temp015);
        }

        private static void PopulateDeterminant_L()
        {
            var Temp020 = new List<object>();
            Temp020.Add(new Matrix(Temp021));
            Temp020.Add(Temp022);
            Determinant_L.Add(Temp020);
            var Temp023 = new List<object>();
            Temp023.Add(new Matrix(Temp024));
            Temp023.Add(Temp025);
            Determinant_L.Add(Temp023);
            var Temp026 = new List<object>();
            Temp026.Add(new Matrix(Temp027));
            Temp026.Add(Temp028);
            Determinant_L.Add(Temp026);
            var Temp029 = new List<object>();
            Temp029.Add(new Matrix(Temp030));
            Temp029.Add(Temp031);
            Determinant_L.Add(Temp029);
        }

        private static void PopulateGramSchmidt_L()
        {
            var Temp032 = new List<object>();
            Temp032.Add(new Matrix(Temp033));
            Temp032.Add(new Matrix(Temp034));
            Temp032.Add(new Matrix(Temp035));
            GramSchmidt_L.Add(Temp032);
            var Temp036 = new List<object>();
            Temp036.Add(new Matrix(Temp037));
            Temp036.Add(new Matrix(Temp038));
            Temp036.Add(new Matrix(Temp039));
            GramSchmidt_L.Add(Temp036);
            var Temp040 = new List<object>();
            Temp040.Add(new Matrix(Temp041));
            Temp040.Add(new Matrix(Temp042));
            Temp040.Add(new Matrix(Temp043));
            GramSchmidt_L.Add(Temp040);
            var Temp044 = new List<object>();
            Temp044.Add(new Matrix(Temp045));
            Temp044.Add(new Matrix(Temp046));
            Temp044.Add(new Matrix(Temp047));
            GramSchmidt_L.Add(Temp044);
        }

        private static void InitAllLists()
        {
            PopulateSquareSubMatrix_L();
            PopulateDeterminant_L();
            PopulateGramSchmidt_L();
        }

        // result lists for the tests
        private static readonly List<bool> res_SquareSubMatrix_L = new List<bool>();
        private static readonly List<bool> res_Determinant_L = new List<bool>();
        private static readonly List<bool> res_GramSchmidt_L = new List<bool>();

        //Routine that runs the differents tests, to be called by Main
        private static void RunTests()
        {
            foreach (List<object> item in SquareSubMatrix_L)
            {
                res_SquareSubMatrix_L.Add(TestSquareSubMatrix((Matrix)item[0], (int)item[1], (int)item[2], (Matrix)item[3]));
            }

            foreach (List<object> item in Determinant_L)
            {
                res_Determinant_L.Add(TestDeterminant((Matrix)item[0], (double)item[1]));
            }

            foreach (List<object> item in GramSchmidt_L)
            {
                res_GramSchmidt_L.Add(TestGramSchmidt((Matrix)item[0], (Matrix)item[1], (Matrix)item[2]));
            }

        }

        // call the reporting routine for each case tested
        private static void PrintReport()
        {
            PrintSummary(res_SquareSubMatrix_L, "Matrix.SquareSubMatrix");
            PrintSummary(res_Determinant_L, "Matrix.Determinant");
            PrintSummary(res_GramSchmidt_L, "Matrix.GramSchmidt");
        }

    }
}

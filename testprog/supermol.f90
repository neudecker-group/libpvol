module supermol
!> large geometry for advanced testing of gradient
  use iso_fortran_env,only:wp => real64

!&<
  integer,parameter :: testnat = 226
  !> Atomtypes
  integer, parameter :: testat(testnat) = [6,6,6,6,6,6,6,6,6,6,   &
  &  6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6, & 
  &  6,6,6,6,6,6,6,6,6,6,6,1,15,1,16,1,1,17,7,1,16,9,1,1,9,1,17,  &
  &  1,7,1,1,7,15,1,1,15,9,1,1,17,1,6,5,1,1,1,1,1,17,1,5,15,1,15, &
  &  1,1,1,1,16,6,9,1,1,8,1,1,16,1,8,16,1,1,15,3,3,1,1,1,1,8,     &
  &  1,17,1,1,8,8,1,6,1,17,1,1,1,1,5,6,5,1,1,1,1,7,1,1,1,1,6,1,3, &
  &  17,1,3,7,9,1,1,1,1,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,6,1,  &
  &  1,1,1,1,1,1,1,1,1,1,6,1,1,1,1,1,6,1,1,1,1,1,1,1,1,1,1,1,1,   &
  &  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

  !> geometry in bohr
real(wp),parameter :: testxyz(3,testnat) =    reshape(&
  [ -9.9390583161_wp,	6.1477393549_wp,	5.3760069171_wp,  &
&   -7.5065206674_wp,	4.4465235226_wp,	5.5738233749_wp,  &
&   -4.8906158997_wp,	5.1596368902_wp,	4.1564126027_wp,  &
&   -3.2679771782_wp,	2.5389168726_wp,	3.7439579532_wp,  &
&   -0.370387564_wp,	2.6447830441_wp,	2.8757760089_wp,  &
&    1.291548969_wp, 	0.0764126946_wp, 	2.1265149262_wp,  &
&    3.947450856_wp,	1.2901951227_wp,	1.2361772787_wp,  &
&    6.2701133085_wp,	0.2116766175_wp,	-0.265667071_wp,  &
&    8.5592533781_wp,	2.2868642024_wp,	-0.5008048796_wp, &
&    11.3328603663_wp, 1.2569982145_wp, 0.6946935397_wp, &
&    13.4456252615_wp, 3.50199821_wp,	  1.3558705524_wp, &
&   -0.4031320567_wp,	-1.6405245291_wp,	0.2133160806_wp, &
&   -0.3622625972_wp,	-4.6954786106_wp,	-0.1047575174_wp, &
&    1.2647569042_wp,	-5.8564678334_wp,	-2.4331564851_wp, &
&    1.7779547572_wp,	-8.9416158798_wp,	-2.274473118_wp, &
&    0.4763827542_wp,	-11.4131181204_wp,	-3.6168186193_wp, &
&   -2.5061823916_wp,	-11.7203384689_wp,	-3.2939694404_wp, &
&   -4.1326982856_wp,	-13.6155256995_wp,	-4.9631500995_wp, &
&   -5.3863609039_wp,	 6.8736017724_wp,	  1.6258145083_wp, &
&   -5.7459190198_wp,	9.8589601825_wp,	2.0879647809_wp, &
&   -6.2362221928_wp,	11.6205038744_wp,	-0.4025935028_wp, &
&   -3.7311401437_wp,	12.9517550898_wp,	-1.4668351349_wp, &
&   -3.7847442967_wp,	13.5118955373_wp,	-4.3632800151_wp, &
&   -8.7592096951_wp,	13.3794817376_wp,	-0.03772333_wp, &
&   -9.5025374288_wp,	15.2495122107_wp,	-2.3550996049_wp, &
&    8.2048095408_wp,	4.2399385157_wp,	-2.9356846994_wp, &
&    9.4339414162_wp,	7.0185955565_wp,	-2.5212356437_wp, &
&    9.3037010999_wp,	2.8586739474_wp,	-5.2413525745_wp, &
&    8.0695240833_wp,	8.9648774036_wp, 	-0.6348407636_wp, &
&    9.8934255542_wp,	8.3887138914_wp,	-5.123672597_wp, &
&    15.7777899507_wp, 	4.0028667489_wp, 	-0.4780754871_wp, &
&    14.7515191076_wp,	3.352716207_wp,	4.0859549806_wp, &
&    12.7800567016_wp, -1.2408674768_wp, 	-0.4654692456_wp, &
&   -6.8765284001_wp,	-13.9846268545_wp,	-3.9417572211_wp, &
&   -8.2318979897_wp,	-11.5488150021_wp,	-3.2863378804_wp, &
&   -3.0964633453_wp,	-5.8823590341_wp,	0.1872208679_wp, &
&   -4.0754523295_wp,	-6.0706175297_wp,	2.9525734928_wp, &
&   -5.1270580628_wp,	-4.6332963197_wp,	-1.5523233411_wp, &
&    11.6461066248_wp,	-4.0222267321_wp,	-0.1701949331_wp, &
&    13.5320670349_wp,	-6.17657664_wp,	-1.0230191185_wp, &
&    12.4916351144_wp,	-8.918907833_wp,	-1.1208256826_wp, &
&    13.5866319333_wp,	-0.9865478177_wp,	-3.2516290166_wp, &
&   -12.4794678765_wp,	4.5957567894_wp,	5.5964435832_wp, &
&   -12.4050492103_wp, 	15.7264219959_wp,	-2.8089178763_wp, &
&   -8.3294965749_wp,	18.017929927_wp,	-2.4353076783_wp, &
&    1.8570522875_wp,	-13.934928209_wp,	-2.5914327727_wp, &
&    1.5090939539_wp,	-15.1043540836_wp,	0.1788442863_wp, &
&    3.2027643646_wp,	-14.09875929_wp,	2.4678761162_wp, &
&    2.5475259901_wp,	-11.421099901_wp,	3.4521612897_wp, &
&    6.1066991298_wp,	-14.4513320728_wp,	2.1228591445_wp, &
&    1.8900567066_wp,	-18.0443330467_wp,	0.0212083697_wp, &
&   -10.0414617617_wp,	6.8164249537_wp,	3.4186402291_wp, &
&   -9.8011203248_wp, 	9.0046389646_wp,	7.547414495_wp,  &
&   -8.0783398574_wp,	2.6852711145_wp,	4.5838968918_wp, &
&   -6.9780646058_wp,	3.7496982464_wp,	8.9631901488_wp, &
&   -3.7614727612_wp,	6.3336292334_wp,	5.4533777125_wp, &
&   -4.392049674_wp,	1.4918128152_wp,	2.3286721171_wp, &
&   -3.1992844488_wp,	0.542163779_wp,	6.4634341352_wp, &
&    1.0617047051_wp,	4.5487622474_wp,	4.3742785287_wp, &
&   -0.4311767174_wp,	3.465379945_wp,	1.0243043941_wp, &
&    1.9176055048_wp,	-1.8873888719_wp,	4.9365252112_wp, &
&    5.0077066981_wp,	2.2308889766_wp,	3.4403040367_wp, &
&    3.402400685_wp,	2.8624976612_wp,	0.0794020371_wp, &
&    5.7654928378_wp,	-0.4480798672_wp,	-2.171596289_wp, &
&    6.9779153748_wp,	-1.9054297664_wp,	1.026495823_wp,  &
&    8.1463769149_wp,	3.5797553059_wp,	1.0123428683_wp, &
&    10.2917425499_wp,	0.3469881432_wp,	3.8607108095_wp, &
&    12.3596894089_wp,	5.248518496_wp,	1.6352715979_wp, &
&   -0.7817337483_wp,	-0.461212573_wp,	-2.2452950573_wp, &
&   -2.158735977_wp,	-1.5090505069_wp,	1.2550389437_wp, &
&    0.6177793156_wp,	-5.4452699356_wp,	1.582772321_wp, &
&    3.6958054542_wp,	-4.5727338463_wp,	-2.0871477884_wp, &
&    0.2062778555_wp,	-4.7765416454_wp,	-5.6834632017_wp, &
&    1.3799619139_wp,	-9.265136976_wp,	-0.2963279967_wp, &
&    4.5164180338_wp,	-9.257104948_wp,	-2.5866770263_wp, &
&    1.2391656767_wp,	-11.2867136075_wp,	-7.1348818298_wp, &
&   -3.5792921231_wp,	-9.4405736952_wp,	-3.8829613405_wp, &
&   -2.862227912_wp,	-12.1304730486_wp,	-1.2741193549_wp, &
&   -4.2395505151_wp,	-12.786850679_wp,	-6.8896673396_wp, &
&   -2.950744563_wp,	-16.7375056568_wp,	-5.4257044292_wp, &
&   -7.1287017774_wp,	6.1389460599_wp,	0.7151850339_wp, &
&   -3.5161568317_wp,	6.4394959716_wp,	-0.4945524785_wp, &
&   -3.7501646405_wp,	11.0702893663_wp,	4.0579037084_wp, &
&   -7.4510382817_wp,	9.9546019526_wp,	3.2147211679_wp, &
&   -6.8158710989_wp,	10.2976445415_wp,	-1.9172966533_wp, &
&   -1.4243468024_wp,	11.4803361075_wp,	-0.9977493988_wp, &
&   -3.4984975152_wp,	14.6742871537_wp,	-0.3244321609_wp, &
&   -3.9306810696_wp,	11.7061705911_wp,	-5.4227063973_wp, &
&   -0.9836876828_wp,	15.107415397_wp,	-5.372308444_wp, &
&   -5.3292029752_wp, 	14.685737705_wp,	-5.0437461735_wp, &
&   -8.5352776234_wp,	15.0703817409_wp,	2.5330462384_wp, &
&   -11.3577490747_wp, 10.8887876348_wp, 	0.3527307032_wp, &
&   -9.0270593167_wp,	14.1786490669_wp,	-4.0800951137_wp, &
&    4.8648654295_wp,	4.9622420572_wp,	-4.1269602592_wp, &
&    11.2726275234_wp,	6.7072065628_wp,	-1.7396584935_wp, &
&    8.4748080888_wp,	0.9421442281_wp,	-5.4269419131_wp, &
&    11.3494053617_wp,	2.7842501809_wp,	-5.1307750481_wp, &
&    8.8772448046_wp,	3.8158670829_wp,	-7.0582335769_wp, &
&    10.1857643607_wp,	11.5947328141_wp,	0.1848176821_wp, &
&    7.1248397169_wp,	7.9896188166_wp,	1.9253944185_wp, &
&    5.9968386678_wp,	10.0185787935_wp,	-1.8345361631_wp, &
&    8.0507364671_wp,	8.8084196563_wp,	-6.0403753953_wp, &
&    11.0277184845_wp,	7.1858098365_wp,	-6.4023512423_wp, &
&    11.2679934286_wp,	10.6318033523_wp,	-4.8570096452_wp, &
&    16.913099163_wp,	2.2516999007_wp,	-0.6526496889_wp, &
&    17.064952095_wp,	5.403905747_wp,	0.3894199946_wp, &
&    15.4102630883_wp,	5.4402319726_wp,	-3.5924924551_wp, &
&    13.371518847_wp,	3.2507088847_wp,	5.6604591486_wp, &
&    16.137592881_wp,	5.5564502419_wp,	4.6223775219_wp, &
&    16.9004870656_wp,	0.6753095662_wp,	4.503636682_wp, &
&    14.5457040571_wp,	-1.3296478322_wp,	0.6279542252_wp, &
&   -6.7828545902_wp,	-15.0934636208_wp, 	-2.1614935392_wp, &
&   -8.8837114183_wp, -15.7829495653_wp, 	-6.2538816982_wp, &
&   -6.3237869767_wp, 	-10.0897670808_wp, 	-0.2720169445_wp, &
&   -11.7840264817_wp, 	-12.3388514721_wp, 	-1.9494115423_wp, &
&   -8.3640463849_wp, 	-10.2707264804_wp, 	-4.9448732689_wp, &
&   -2.8188959521_wp, 	-7.8790106884_wp, 	-0.2017028864_wp, &
&   -2.5445425035_wp, 	-6.7288777375_wp, 	4.2359806829_wp,  &
&   -5.5514676266_wp, 	-7.5768335921_wp, 	2.9341746328_wp,  &
&   -5.1387963088_wp, 	-3.7897682811_wp, 	3.8219647749_wp,  &
&   -5.3525062459_wp, 	-2.6031141598_wp, 	-1.0949902074_wp, &
&   -8.2238771068_wp, 	-5.9259788776_wp, 	-1.2325220637_wp, &
&   -4.6093184871_wp, 	-4.7410176781_wp, 	-3.5572002681_wp, &
&    9.9497056646_wp, 	-4.1589650195_wp, 	-1.3941994958_wp, &
&    10.9021171744_wp, 	-4.5261182377_wp, 	2.3242267334_wp,  &
&    15.797748558_wp, 	-6.0541008693_wp, 	0.3668428555_wp,  &
&    14.0696991029_wp, 	-5.9029600903_wp, 	-3.0166264604_wp, &
&    11.3682156459_wp, 	-9.9549213825_wp, 	1.3339296285_wp,  &
&    14.0931941071_wp, 	-10.1666514587_wp, 	-1.6651019303_wp, &
&    10.2364148733_wp, 	-9.1917345918_wp, 	-3.6030774666_wp, &
&    11.9750259163_wp,	-1.3448106463_wp, 	-4.5306244472_wp, &
&    15.1204372185_wp, 	-2.2971405513_wp, 	-3.7636464335_wp, &
&    14.4802445401_wp, 	0.7760769489_wp,  	-3.7116569977_wp, &
&   -12.405597587_wp, 	3.0532105345_wp,  	4.1709165336_wp,  &
&   -14.8876481925_wp, 	6.2655572294_wp,  	4.7860164916_wp,  &
&   -13.0261070743_wp, 	3.3427426648_wp,  	8.1582855487_wp,  &
&   -13.8708486065_wp, 	17.3465404136_wp, 	-0.6908615143_wp, &
&   -13.4448124248_wp, 	13.9488518391_wp, 	-3.1452052123_wp, &
&   -12.6202652471_wp,	16.7629419863_wp, 	-4.6231692018_wp, &
&   -9.1331361117_wp, 	19.2126502672_wp, 	-0.9042835678_wp, &
&   -8.9462307068_wp, 	18.8976912524_wp, 	-4.2451867244_wp, &
&   -5.5907393699_wp, 	18.2558394653_wp, 	-2.3237321932_wp, &
&    3.8942235585_wp, 	-13.9657986248_wp, 	-3.0565166623_wp, &
&    1.2018437592_wp, 	-15.3933002461_wp, 	-3.8851166895_wp, &
&   -0.4911029756_wp, 	-14.877661994_wp, 	0.7563846906_wp,  &
&    2.6387469335_wp, 	-15.2645201268_wp, 	4.1193523117_wp,  &
&    3.3322861053_wp, 	-10.9596181455_wp, 	6.2020832403_wp,  &
&    0.4588501896_wp, 	-11.1847525657_wp, 	3.4196202423_wp,  &
&    4.4571487698_wp, 	-8.5364676934_wp, 	1.7525280616_wp,  &
&    7.1798685537_wp, 	-17.6007440881_wp, 	2.6849624519_wp,  &
&    7.1544165455_wp, 	-13.3565005414_wp, 	3.5569616874_wp,  &
&    7.6382665167_wp, 	-13.4109145087_wp, 	-1.2815074228_wp, &
&    0.0090341425_wp, 	-19.2546750391_wp, 	-1.5789040268_wp, &
&    1.7200750416_wp, 	-19.1545066583_wp, 	2.3914127964_wp,  &
&    3.7782151571_wp, 	-18.4534794905_wp, 	-0.7887749958_wp, &
&   -10.2542351326_wp, 	7.7585212012_wp,   	9.9395736353_wp,  &
&   -12.3325087768_wp, 	9.9490090014_wp,  	7.125126656_wp,   &
&   -5.6041792271_wp, 	5.9731942437_wp,  	9.448934787_wp,   &
&    1.2342965903_wp, 	4.2272563364_wp,  	7.118974183_wp,   &
&    0.3872759796_wp, 	6.3698868883_wp,  	4.0091148053_wp,  &
&    3.7586003478_wp, 	-3.4217240216_wp, 	3.8125343202_wp,  &
&    0.966765643_wp,  	-0.2102245432_wp, 	-3.1040378575_wp, &
&   -1.6914191005_wp,  	1.2576871328_wp,  	-2.1365696954_wp, &
&    4.5629080239_wp, 	-5.1892103186_wp, 	-0.4269192892_wp, &
&    4.9289233085_wp, 	-4.8280150064_wp, 	-3.6120162146_wp, &
&   -1.9358281211_wp, 	-6.3296827346_wp, 	-6.1995202451_wp, &
&    2.1452096036_wp, 	-6.0907006221_wp, 	-7.0988880682_wp, &
&    5.0898543967_wp, 	-9.4054067489_wp, 	-4.4661062189_wp, &
&    5.3012714656_wp, 	-10.5928368714_wp, 	-1.3889692326_wp, &
&   -0.8592552859_wp, 	-9.7814069097_wp, 	-8.0070571201_wp, &
&    0.1799865646_wp, 	-13.6689611648_wp, 	-7.9428160081_wp, &
&   -1.5643217024_wp,  	11.1683175101_wp, 	3.6519011273_wp,  &
&   -4.4635298118_wp, 	11.7685563273_wp,  	6.0546707806_wp,  &
&   -1.0035850865_wp, 	10.2671789384_wp, 	-2.4985474314_wp, &
&    0.1261117179_wp, 	12.69542119_wp,    	-0.8305772435_wp, &
&   -6.8727490778_wp,  	16.5289727131_wp, 	2.7037957819_wp,  &
&   -10.3696173305_wp, 	15.0061039823_wp,  	4.919252573_wp,   &
&   -11.4991089895_wp, 	10.1647216563_wp, 	-2.2917147766_wp, &
&   -13.6237955526_wp, 	12.3643818738_wp, 	0.681447305_wp,   &
&    3.6457307088_wp, 	6.1145454691_wp,  	-1.9784524616_wp, &
&    3.8479145346_wp, 	2.4561527933_wp,  	-4.2782361017_wp, &
&    11.9872966714_wp, 	10.0500695569_wp, 	1.3888728936_wp,  &
&    11.3097609187_wp, 	11.4369374659_wp, 	-6.5500115388_wp, &
&    17.6906074072_wp, 	4.4193043549_wp,  	-4.5037018525_wp, &
&    14.901812859_wp, 	6.9406625947_wp,  	4.9035602722_wp,  &
&    15.1727026372_wp, 	-0.7756974488_wp, 	5.9087906777_wp,  &
&   -9.1470344292_wp, 	-13.7210373057_wp, 	-8.0467550218_wp, &
&   -11.2367302897_wp, 	-15.4823556399_wp, 	-4.8708048233_wp, &
&   -6.6796823216_wp, 	-4.2896426393_wp, 	5.9483640943_wp,  &
&    12.4011098897_wp, 	-4.3213119132_wp, 	3.4291574255_wp,  &
&    15.5112526174_wp, 	-6.8051041767_wp, 	2.0539986719_wp,  &
&   -14.9605898496_wp, 	7.2167097207_wp,  	2.7709905343_wp,  &
&   -16.618757072_wp,  	6.5406702449_wp,  	6.1693497038_wp,  &
&   -14.0501181067_wp, 	16.6028086228_wp, 	1.3968148163_wp,  &
&   -15.1233619572_wp, 	19.9937959698_wp, 	-1.3589343215_wp, &
&   -5.0322281849_wp, 	18.4140893069_wp, 	-0.4426505034_wp, &
&   -5.0791668125_wp, 	19.9835303882_wp, 	-3.1349515314_wp, &
&   -1.7644295172_wp, 	-19.0747421221_wp,	-0.7246381566_wp, &
&    0.4223870115_wp, 	-21.182477177_wp, 	-1.7023097522_wp, &
&   -1.6026149437_wp, 	7.034086301_wp,   	0.0759062537_wp,  &
&   -4.052313024_wp,  	7.4106858415_wp,  	-2.2592670458_wp, &
&   -3.6473106684_wp, 	4.4548785696_wp,  	-1.0602361965_wp, &
&    8.6857887987_wp, 	7.1977088364_wp,  	3.0781026431_wp,  &
&    5.5866732974_wp, 	6.6050279082_wp,  	1.7031388343_wp,  &
&    6.2808139268_wp, 	9.5712005009_wp,  	3.0191953115_wp,  &
&    12.7541607784_wp, 	-9.8858349855_wp, 	2.9017176459_wp,  &
&    10.8490965335_wp, 	-11.9639875448_wp, 	1.0447375182_wp,  &
&    9.6299187436_wp, 	-8.9172271897_wp,  	1.8732855979_wp,  &
&   -11.5743316141_wp, 	1.9208477844_wp,  	8.634065523_wp,   &
&   -14.8369451657_wp, 	2.2806429698_wp,  	8.059601733_wp,   &
&   -13.1964771177_wp, 	4.7504116751_wp,  	9.6982151608_wp,  &
&    2.4109769805_wp, 	-12.3451650768_wp, 	7.4795548146_wp,  &
&    2.7214486758_wp, 	-9.0326626574_wp, 	6.7609416235_wp,  &
&    5.4050858053_wp, 	-11.0555780232_wp, 	6.4554124967_wp,  &
&    1.9678723615_wp, 	2.3473360038_wp, 	  7.63837929_wp,    &
&   -0.606488007_wp,  	4.5882618041_wp,  	8.0483041907_wp,  &
&    2.5859261488_wp,  	5.6416377193_wp,  	7.8741410353_wp,  &
&   -11.9417827276_wp, 	13.6545014365_wp,  	4.6985301012_wp,  &
&   -11.1846358255_wp, 	16.9170214826_wp,  	5.1915541898_wp,  &
&   -9.2556677968_wp, 	14.505107747_wp,  	6.6196172219_wp,  &
&   -7.5442961338_wp,  	-2.4913690425_wp,  	6.5941240335_wp,  &
&   -5.5386683293_wp,  	-5.05046818_wp,   	7.5444036179_wp,  &
&   -8.2799435713_wp,  	-5.5654316813_wp,  	5.4646533424_wp,  &
&   -14.817153546_wp,  	20.4732258699_wp,  	-3.3765196849_wp, &
&   -14.2629762937_wp, 	21.4842800624_wp,  	-0.1632507802_wp, &
&   -17.182451673_wp, 	19.9000447109_wp, 	-0.9803347078_wp], &
&  shape(testxyz))

  !> testpressure in GPa
  real(wp), parameter :: testpressure = 10.0_wp

  !> proberad in GPa
  real(wp), parameter :: testproberad = 0.0_wp

  ! Volume using unscaled Bondi Radii
  !> to large to compute Volume numerically
  !public :: testnat
  !public :: testat
  !public :: testxyz
  !public :: testpressure
  !public :: testproberad
  !public :: testgrad
  !public :: writetestcoord

  character(len=2),parameter :: testelem(18) = ['h ','he','li','be', &
      &                        'b ','c ','n ','o ', 'f ', 'ne', 'na', 'mg', 'al', 'si', 'p ', 's ', 'cl', 'ar']

contains

  subroutine writetestcoord

    integer :: ich,i

    open (newunit=ich,file='coord')
    write (ich,'(a)') '$coord'
    do i = 1,testnat
      write (ich,'(3f20.14,2x,a2)') testxyz(1:3,i),testelem(testat(i))
    end do
    write (ich,'(a)') '$end'
    close (ich)

  end subroutine writetestcoord

end module supermol

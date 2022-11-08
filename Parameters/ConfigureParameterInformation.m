[pNames, param] = IQMparameters(model);
AmountParameterOpti = 73;

[row, column] = size(param);
ParametersInformation(row+1:end,:) = [];

ParametersInformation(1:row,1) = strcat(pNames);
ParametersInformation(1:row,4) = array2table(param);

lb = 0.75;
ub = 1.5;

ParametersInformation((1:AmountParameterOpti),2) = array2table(param(1:AmountParameterOpti)*lb);
ParametersInformation((1:AmountParameterOpti),3) = array2table(param(1:AmountParameterOpti)*ub);
%% Parameters fixed to a specific intervall

% aa ratio
ParametersInformation.LowerBound(1:AmountParameterOpti) = [0.0116527644618584;0.231727017901948;0.00285912408468584;0.0114681139626515;0.0114444108366990;0.00322171180333014;199.585365659147;2.85093161519303;20.1577650000000;1.59874746763680;0.479959660852995;0.530272709027517;0.345206134016510;2.71033082129471;0.0346437224898594;0.00492353833973691;0.667603800224545;2.15911023487725;28.5371568094764;7.78634250000000e-06;1.21791107422657;7.02281870761048e-05;0.000171623429601669;0.0334091962500000;1.74226912500000;3.66388875000000;0.913803193971135;0.000486628896298269;0.00373166922180501;0.407584714293327;15.0231612074152;0.291656164619505;127.891636845264;40.5100425000000;1.79591751462143;0.937059143153920;0.0120955521153968;11.6379917932246;13.8249937500000;0.434683099595435;0.347500603269185;0.599325359683343;0.00582143954437330;0.429085312500000;0.101943666811664;0.111583819250565;0.0913560375000000;0.0277982357156061;0.000416249027328152;2.90476987500000e-05;0.000519641836758583;0.676784512500000;0.00417817441230182;0.465533250000000;0.0173705045142078;3.27789886630369;132.505386328950;3.80005752605088;258.631164884812;16.1770748451646;0.0336766687500000;0.0287459209990777;0.708861003456152;24.4315895954794;4.98206070532111e-05;0.00709865887500000;0.861422428803124;0.0321477150000000;1728.75355295055;392.832285081669;0.0946872599608725;0.0694986138170292;0.0631950783665003];
ParametersInformation.UpperBound(1:AmountParameterOpti) = [0.0184682310863370;0.512035650000000;0.00630359100000000;0.0193873574156005;0.0244756461933206;0.00555009005979407;405.558898612496;5.48283023361430;30.5240870692219;3.20427324185572;0.793543275000000;0.703166284746729;0.481394615240764;3.89672405393854;0.0765497250000000;0.00793549758326067;0.780202328147862;5.25000000000000;60.7249125000000;1.69709514394991e-05;2.39610349468232;0.000122170473335141;0.000358156575000000;0.0738485344496604;2.18332802633878;5.64939362623991;1.64294550000000;0.00101219037580826;0.00535253881673252;0.714327075000000;33.1652338769216;0.517256733318133;277.578215934311;89.5485150000000;3.96596659114417;1.30763595533976;0.0236550296927188;25.5723387342206;30.5605125000000;0.801649177084107;0.708679125000000;1.13673841788036;0.0128584575000000;0.948487241032062;0.944960912960254;0.213520237547503;0.195214822033223;0.0572948370606018;0.000886529700000000;6.40559346342997e-05;0.000792086128266826;1.49385496329594;0.00575710396198530;0.713014210974871;0.0233717898001979;5.27348975684954;224.652753713648;4.41163885787888;565.449312394239;35.6136971947612;0.0634772770650199;0.0549231039674964;1.40317193700442;51.2718612651326;8.02362463360413e-05;0.0139780585877922;1.88886600000000;0.0602734727245523;2354.46153840959;524.621715840200;0.188969223898596;0.115499992636239;0.104802632375319];

/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"lexer.h"

void lexer::ctrlsend()
{
    int n;
    int ii,dd;
    
    ii=dd=0;

    for(n=0;n<ctrlsize;++n)
    {
    ictrl[n]=0;
    dctrl[n]=0.0;
    }
    
    ictrl[ii] = A10;
	ii++;
	ictrl[ii] = A209;
	ii++;
    ictrl[ii] = A210;
	ii++;
	ictrl[ii] = A211;
	ii++;
	ictrl[ii] = A212;
	ii++;
    ictrl[ii] = A214;
	ii++;
    ictrl[ii] = A215;
	ii++;
    ictrl[ii] = A216;
	ii++;
    ictrl[ii] = A219;
	ii++;
	ictrl[ii] = A220;
	ii++;
	ictrl[ii] = A221;
	ii++;
    ictrl[ii] = A222;
	ii++;
    dctrl[dd] = A223;
	dd++;
	ictrl[ii] = A230;
	ii++;
    ictrl[ii] = A240;
	ii++;
	ictrl[ii] = A241;
	ii++;
	ictrl[ii] = A242;
	ii++;
    ictrl[ii] = A243;
	ii++;
    ictrl[ii] = A244;
	ii++;
    dctrl[dd] = A244_val;
	dd++;
    ictrl[ii] = A245;
	ii++;
    dctrl[dd] = A245_val;
	dd++;
    ictrl[ii] = A246;
	ii++;
    dctrl[dd] = A247;
	dd++;
    ictrl[ii] = A248;
	ii++;
    dctrl[dd] = A249;
	dd++;
    ictrl[ii] = A260;
	ii++;
    ictrl[ii] = A300;
	ii++;
    ictrl[ii] = A310;
	ii++;
    ictrl[ii] = A311;
	ii++;
    ictrl[ii] = A312;
	ii++;
    ictrl[ii] = A313;
	ii++;
    ictrl[ii] = A320;
	ii++;
    ictrl[ii] = A321;
	ii++;
    ictrl[ii] = A410;
	ii++;
    dctrl[dd] = A440;
	dd++;


    ictrl[ii] = M10;
    ii++;

    ictrl[ii] = I10;
    ii++;
    ictrl[ii] = I11;
    ii++;
    ictrl[ii] = I12;
    ii++;
    ictrl[ii] = I13;
    ii++;
    ictrl[ii] = I20;
    ii++;
    ictrl[ii] = I21;
    ii++;
	ictrl[ii] = I30;
    ii++;
	ictrl[ii] = I40;
    ii++;
	ictrl[ii] = I41;
    ii++;
    dctrl[dd] = I55;
    dd++;
    ictrl[ii] = I56;
    ii++;
    dctrl[dd] = I58_1;
    dd++;
    dctrl[dd] = I58_2;
    dd++;
    ictrl[ii] = I230;
    ii++;
    dctrl[dd] = I231;
    dd++;
    dctrl[dd] = I232;
    dd++;
	
    ictrl[ii] = B10;
    ii++;
    ictrl[ii] = B19;
    ii++;
    ictrl[ii] = B20;
    ii++;
	ictrl[ii] = B26;
    ii++;
	ictrl[ii] = B28;
    ii++;
	dctrl[dd] = B29;
    dd++;
    ictrl[ii] = B30;
    ii++;
    dctrl[dd] = B50;
    dd++;
	dctrl[dd] = B51;
    dd++;
	dctrl[dd] = B52;
    dd++;
	dctrl[dd] = B53;
    dd++;
	dctrl[dd] = B54;
    dd++;
	dctrl[dd] = B55;
    dd++;
	dctrl[dd] = B56;
    dd++;
    ictrl[ii] = B60;
    ii++;
    ictrl[ii] = B61;
    ii++;
    ictrl[ii] = B62;
    ii++;
    ictrl[ii] = B63;
    ii++;
    ictrl[ii] = B64;
    ii++;
    dctrl[dd] = B65;
    dd++;
	dctrl[dd] = B66_1;
    dd++;
	dctrl[dd] = B66_2;
    dd++;
    ictrl[ii] = B67;
    ii++;
	ictrl[ii] = B68;
    ii++;
	ictrl[ii] = B69;
    ii++;
	ictrl[ii] = B70;
    ii++;
	ictrl[ii] = B71;
    ii++;
    ictrl[ii] = B74;
    ii++;
	ictrl[ii] = B75;
    ii++;
	ictrl[ii] = B76;
    ii++;
	ictrl[ii] = B77;
    ii++;
	ictrl[ii] = B78;
    ii++;
	dctrl[dd] = B79;
    dd++;
	dctrl[dd] = B80;
    dd++;
	ictrl[ii] = B81;
    ii++;
	dctrl[dd] = B81_1;
    dd++;
	dctrl[dd] = B81_2;
    dd++;
    ictrl[ii] = B82;
    ii++;
	dctrl[dd] = B83;
    dd++;
	ictrl[ii] = B84;
    ii++;
    ictrl[ii] = B85;
    ii++;
	ictrl[ii] = B86;
    ii++;
	ictrl[ii] = B87;
    ii++;
	dctrl[dd] = B87_1;
    dd++;
	dctrl[dd] = B87_2;
    dd++;
	dctrl[dd] = B88;
    dd++;
    ictrl[ii] = B89;
    ii++;
	ictrl[ii] = B90;
    ii++;
    ictrl[ii] = B91;
    ii++;
    dctrl[dd] = B91_1;
    dd++;
    dctrl[dd] = B91_2;
    dd++;
    dctrl[dd] = B91_3;
    dd++;
    ictrl[ii] = B92;
    ii++;
    ictrl[ii] = B93;
    ii++;
    dctrl[dd] = B93_1;
    dd++;
    dctrl[dd] = B93_2;
    dd++;
    dctrl[dd] = B93_3;
    dd++;
    dctrl[dd] = B96_1;
    dd++;
    dctrl[dd] = B96_2;
    dd++;
    dctrl[dd] = B96_3;
    dd++;
    dctrl[dd] = B97;
    dd++;
    ictrl[ii] = B98;
    ii++;
    ictrl[ii] = B99;
    ii++;
    ictrl[ii] = B101;
    ii++;
    dctrl[dd] = B102;
    dd++;
	dctrl[dd] = B103;
    dd++;
	dctrl[dd] = B104;
    dd++;
	dctrl[dd] = B105_1;
    dd++;
	dctrl[dd] = B105_2;
    dd++;
	dctrl[dd] = B105_3;
    dd++;
	ictrl[ii] = B106;
    ii++;
	ictrl[ii] = B107;
    ii++;
    ictrl[ii] = B108;
    ii++;
	ictrl[ii] = B109;
    ii++;
	dctrl[dd] = B110_d;
    dd++;
	ictrl[ii] = B110;
    ii++;
	dctrl[dd] = B111_zs;
    dd++;
	dctrl[dd] = B111_ze;
    dd++;
    dctrl[dd] = B112_zs;
    dd++;
    dctrl[dd] = B112_z2;
    dd++;
    dctrl[dd] = B112_ze;
    dd++;
    ictrl[ii] = B115;
    ii++;
    ictrl[ii] = B116;
    ii++;
    dctrl[dd] = B117;
    dd++;
	dctrl[dd] = B118;
    dd++;
    dctrl[dd] = B119;
    dd++;
    dctrl[dd] = B120;
    dd++;
	ictrl[ii] = B121;
    ii++;
    dctrl[dd] = B122;
    dd++;
    dctrl[dd] = B123;
    dd++;
	dctrl[dd] = B126;
    dd++;
    ictrl[ii] = B130;
    ii++;
    dctrl[dd] = B131;
    dd++;
    dctrl[dd] = B132_s;
    dd++;
    dctrl[dd] = B132_e;
    dd++;
    ictrl[ii] = B133;
    ii++;
    dctrl[dd] = B134;
    dd++;
    dctrl[dd] = B135;
    dd++;
    ictrl[ii] = B136;
    ii++;
    dctrl[dd] = B140_1;
    dd++;
    dctrl[dd] = B140_2;
    dd++;
    dctrl[dd] = B140_3;
    dd++;
    ictrl[ii] = B160;
    ii++;
    ictrl[ii] = B180;
    ii++;
    dctrl[dd] = B181_1;
    dd++;
    dctrl[dd] = B181_2;
    dd++;
    dctrl[dd] = B181_3;
    dd++;
    dctrl[dd] = B182_1;
    dd++;
    dctrl[dd] = B182_2;
    dd++;
    dctrl[dd] = B182_3;
    dd++;
    dctrl[dd] = B183_1;
    dd++;
    dctrl[dd] = B183_2;
    dd++;
    dctrl[dd] = B183_3;
    dd++;
	ictrl[ii] = B191;
    ii++;
    dctrl[dd] = B191_1;
    dd++;
    dctrl[dd] = B191_2;
    dd++;
	dctrl[dd] = B191_3;
    dd++;
	dctrl[dd] = B191_4;
    dd++;
	ictrl[ii] = B192;
    ii++;
    dctrl[dd] = B192_1;
    dd++;
    dctrl[dd] = B192_2;
    dd++;
	dctrl[dd] = B192_3;
    dd++;
	dctrl[dd] = B192_4;
    dd++;
	dctrl[dd] = B194_s;
    dd++;
	dctrl[dd] = B194_e;
    dd++;
    ictrl[ii] = B210;
    ii++;
	ictrl[ii] = B240;
    ii++;
	ictrl[ii] = B241;
    ii++;
	ictrl[ii] = B242;
    ii++;
	ictrl[ii] = B243;
    ii++;
    dctrl[dd] = B260;
    dd++;
    dctrl[dd] = B264;
    dd++;
    dctrl[dd] = B267;
    dd++;
    ictrl[ii] = B269;
    ii++;
	ictrl[ii] = B270;
    ii++;
    ictrl[ii] = B274;
    ii++;
    ictrl[ii] = B281;
    ii++;
    ictrl[ii] = B291;
    ii++;
    ictrl[ii] = B295;
    ii++;
    ictrl[ii] = B308;
    ii++;
    dctrl[dd] = B309;
    dd++;
    ictrl[ii] = B310;
    ii++;
    
	
	dctrl[dd] = C1;
    dd++;
    dctrl[dd] = C2;
    dd++;
	dctrl[dd] = C3;
    dd++;
	dctrl[dd] = C4;
    dd++;
	dctrl[dd] = C5;
    dd++;
    ictrl[ii] = C10;
    ii++;
	ictrl[ii] = C15;
    ii++;
	ictrl[ii] = C20;
    ii++;
    dctrl[dd] = C50_1;
    dd++;
    dctrl[dd] = C50_2;
    dd++;
    dctrl[dd] = C51;
    dd++;
    dctrl[dd] = C52;
    dd++;
    dctrl[dd] = C53;
    dd++;
    dctrl[dd] = C54;
    dd++;
    dctrl[dd] = C55;
    dd++;
    dctrl[dd] = C56;
    dd++;
    dctrl[dd] = C57_1;
    dd++;
    dctrl[dd] = C57_2;
    dd++;
    dctrl[dd] = C57_3;
    dd++;
    dctrl[dd] = C57_4;
    dd++;
    dctrl[dd] = C58_1;
    dd++;
    dctrl[dd] = C58_2;
    dd++;
    dctrl[dd] = C58_3;
    dd++;
    dctrl[dd] = C58_4;
    dd++;
	ictrl[ii] = C75;
    ii++;

    ictrl[ii] = D10;
	ii++;
	ictrl[ii] = D11;
	ii++;
    ictrl[ii] = D20;
	ii++;
	ictrl[ii] = D21;
	ii++;
	ictrl[ii] = D22;
	ii++;
	ictrl[ii] = D23;
	ii++;
    ictrl[ii] = D24;
	ii++;
	dctrl[dd] = D29;
	dd++;
    ictrl[ii] = D30;
	ii++;
    ictrl[ii] = D31;
	ii++;
	ictrl[ii] = D32;
	ii++;
	ictrl[ii] = D33;
	ii++;
	ictrl[ii] = D34;
	ii++;
	dctrl[dd] = D35;
	dd++;
    ictrl[ii] = D36;
	ii++;
    ictrl[ii] = D38;
	ii++;
	

    ictrl[ii] = F10;
	ii++;
    ictrl[ii] = F11;
	ii++;
	dctrl[dd] = F19;
	dd++;
    ictrl[ii] = F30;
	ii++;
    ictrl[ii] = F31;
	ii++;
    ictrl[ii] = F32;
	ii++;
    dctrl[dd] = F33;
	dd++;
    ictrl[ii] = F34;
	ii++;
    ictrl[ii] = F35;
	ii++;
	ictrl[ii] = F36;
	ii++;
	dctrl[dd] = F39;
	dd++;
    ictrl[ii] = F40;
	ii++;
    ictrl[ii] = F41;
	ii++;
	dctrl[dd] = F42;
	dd++;
    dctrl[dd] = F43;
	dd++;
    ictrl[ii] = F44;
	ii++;
    dctrl[dd] = F45;
	dd++;
    ictrl[ii] = F46;
	ii++;
    ictrl[ii] = F47;
	ii++;
    ictrl[ii] = F48;
	ii++;
	ictrl[ii] = F49;
	ii++;
    ictrl[ii] = F50;
	ii++;
    ictrl[ii] = F50_flag;
	ii++;
    dctrl[dd] = F51;
	dd++;
    dctrl[dd] = F52;
	dd++;
    dctrl[dd] = F53;
	dd++;
    dctrl[dd] = F54;
	dd++;
    dctrl[dd] = F55;
	dd++;
    dctrl[dd] = F56;
	dd++;
    dctrl[dd] = F57_1;
	dd++;
    dctrl[dd] = F57_2;
	dd++;
    dctrl[dd] = F57_3;
	dd++;
    dctrl[dd] = F57_4;
	dd++;
    dctrl[dd] = F58_1;
	dd++;
    dctrl[dd] = F58_2;
	dd++;
    dctrl[dd] = F58_3;
	dd++;
    dctrl[dd] = F58_4;
	dd++;
    dctrl[dd] = F59_xm;
	dd++;
    dctrl[dd] = F59_ym;
	dd++;
    dctrl[dd] = F59_zs;
	dd++;
    dctrl[dd] = F59_ze;
	dd++;
    dctrl[dd] = F59_r;
	dd++;
    dctrl[dd] = F60;
	dd++;
    dctrl[dd] = F61;
	dd++;
    dctrl[dd] = F62;
	dd++;
    dctrl[dd] = F63;
	dd++;
	ictrl[ii] = F64;
	ii++;
	dctrl[dd] = F64_xs;
	dd++;
	dctrl[dd] = F64_ys;
	dd++;
	dctrl[dd] = F64_zs;
	dd++;
	dctrl[dd] = F64_alpha;
	dd++;
    ictrl[ii] = F70;
	ii++;
	ictrl[ii] = F71;
	ii++;
	ictrl[ii] = F72;
	ii++;
    ictrl[ii] = F80;
	ii++;
    dctrl[dd] = F84;
	dd++;
    ictrl[ii] = F85;
	ii++;
	ictrl[ii] = F101;
	ii++;
    ictrl[ii] = F150;
	ii++;
	ictrl[ii] = F151;
	ii++;
	ictrl[ii] = F300;
	ii++;
	ictrl[ii] = F305;
	ii++;
	ictrl[ii] = F310;
	ii++;
	dctrl[dd] = F321;
	dd++;
	dctrl[dd] = F322;
	dd++;
	dctrl[dd] = F323;
	dd++;
	ictrl[ii] = F350;
	ii++;
	dctrl[dd] = F360;
	dd++;
	dctrl[dd] = F361;
	dd++;
	dctrl[dd] = F362;
	dd++;
    ictrl[ii] = F369;
	ii++;
	ictrl[ii] = F370;
	ii++;
	ictrl[ii] = F371;
	ii++;
	ictrl[ii] = F374;
	ii++;
    ictrl[ii] = F375;
	ii++;
    ictrl[ii] = F378;
	ii++;
    ictrl[ii] = F379;
	ii++;
	dctrl[dd] = F380;
	dd++;
	dctrl[dd] = F381;
	dd++;
	dctrl[dd] = F382;
	dd++;
	ictrl[ii] = F390;
	ii++;
	ictrl[ii] = F391;
	ii++;
    ictrl[ii] = F394;
	ii++;
    ictrl[ii] = F395;
	ii++;
    ictrl[ii] = F398;
	ii++;
    ictrl[ii] = F399;
	ii++;

	
	
    ictrl[ii] = G1;
	ii++;
    ictrl[ii] = G10;
	ii++;
    ictrl[ii] = G11;
	ii++;
    ictrl[ii] = G12;
	ii++;
    ictrl[ii] = G20;
	ii++;
    ictrl[ii] = G21;
	ii++;
    ictrl[ii] = G22;
	ii++;
    ictrl[ii] = G30;
	ii++;
    ictrl[ii] = G40;
	ii++;
    ictrl[ii] = G50;
	ii++;
    ictrl[ii] = G52;
	ii++;
    ictrl[ii] = G60;
	ii++;
    ictrl[ii] = G61;
	ii++;
	ictrl[ii] = G81;
	ii++;
	ictrl[ii] = G95;
	ii++;

    dctrl[dd] = H1;
	dd++;
    dctrl[dd] = H2;
	dd++;
    ictrl[ii] = H10;
	ii++;
    dctrl[dd] = H50_1;
	dd++;
    dctrl[dd] = H50_2;
	dd++;
    dctrl[dd] = H51;
	dd++;
    dctrl[dd] = H52;
	dd++;
    dctrl[dd] = H53;
	dd++;
    dctrl[dd] = H54;
	dd++;
    dctrl[dd] = H55;
	dd++;
    dctrl[dd] = H56;
	dd++;
    dctrl[dd] = H57_1;
	dd++;
    dctrl[dd] = H57_2;
	dd++;
    dctrl[dd] = H57_3;
	dd++;
    dctrl[dd] = H57_4;
	dd++;
    dctrl[dd] = H58_1;
	dd++;
    dctrl[dd] = H58_2;
	dd++;
    dctrl[dd] = H58_3;
	dd++;
    dctrl[dd] = H58_4;
	dd++;

    ictrl[ii] = N5;
    ii++;
    ictrl[ii] = N8;
    ii++;
    ictrl[ii] = N9;
    ii++;
	ictrl[ii] = N10;
    ii++;
    ictrl[ii] = N11;
    ii++;
    ictrl[ii] = N12;
    ii++;
    ictrl[ii] = N13;
    ii++;
    ictrl[ii] = N14;
    ii++;
	ictrl[ii] = N15;
    ii++;
	ictrl[ii] = N16;
    ii++;
	dctrl[dd] = N17;
    dd++;
	dctrl[dd] = N18;
    dd++;
    ictrl[ii] = N40;
    ii++;
    dctrl[dd] = N41;
    dd++;
	ictrl[ii] = N42;
    ii++;
    dctrl[dd] = N43;
    dd++;
	dctrl[dd] = N44;
    dd++;
    ictrl[ii] = N45;
    ii++;
    ictrl[ii] = N46;
    ii++;
    dctrl[dd] = N47;
    dd++;
    ictrl[ii] = N48;
    ii++;
    dctrl[dd] = N49;
    dd++;
    ictrl[ii] = N50;
    ii++;
    dctrl[dd] = N51;
    dd++;
    ictrl[ii] = N52;
    ii++;
    dctrl[dd] = N53;
    dd++;
    dctrl[dd] = N54;
    dd++;
    dctrl[dd] = N55;
    dd++;
    dctrl[dd] = N56;
    dd++;
    ictrl[ii] = N57_1;
    ii++;
    ictrl[ii] = N57_2;
    ii++;
    dctrl[dd] = N57_3;
    dd++;
    ictrl[ii] = N58;
    ii++;
    ictrl[ii] = N60;
    ii++;
    dctrl[dd] = N61;
    dd++;
    ictrl[ii] = N71;
    ii++;
    ictrl[ii] = N72;
    ii++;
    ictrl[ii] = N73;
    ii++;
	
	

    ictrl[ii] = P10;
	ii++;
	ictrl[ii] = P11;
	ii++;
    ictrl[ii] = P12;
	ii++;
	ictrl[ii] = P14;
	ii++;
    ictrl[ii] = P15;
	ii++;
	ictrl[ii] = P17;
	ii++;
	ictrl[ii] = P18;
	ii++;
	ictrl[ii] = P19;
	ii++;
    ictrl[ii] = P20;
	ii++;
	ictrl[ii] = P23;
	ii++;
    ictrl[ii] = P24;
	ii++;
    ictrl[ii] = P25;
	ii++;
	ictrl[ii] = P26;
	ii++;
	ictrl[ii] = P27;
	ii++;
	ictrl[ii] = P28;
	ii++;
	ictrl[ii] = P29;
	ii++;
    dctrl[dd] = P30;
	dd++;
	dctrl[dd] = P34;
	dd++;
    ictrl[ii] = P35;
	ii++;
    ictrl[ii] = P40;
	ii++;
    ictrl[ii] = P41;
	ii++;
	dctrl[dd] = P42;
	dd++;
    ictrl[ii] = P50;
	ii++;
	ictrl[ii] = P51;
	ii++;
    ictrl[ii] = P52;
	ii++;
    ictrl[ii] = P53;
	ii++;
	ictrl[ii] = P54;
	ii++;
	dctrl[dd] = P55;
	dd++;
	ictrl[ii] = P56;
	ii++;
	ictrl[ii] = P59;
	ii++;
	ictrl[ii] = P61;
	ii++;
	ictrl[ii] = P62;
	ii++;
	ictrl[ii] = P66;
	ii++;
	ictrl[ii] = P67;
	ii++;
    ictrl[ii] = P71;
	ii++;
    ictrl[ii] = P75;
	ii++;
	ictrl[ii] = P78;
	ii++;
    ictrl[ii] = P79;
	ii++;
    ictrl[ii] = P81;
	ii++;
    dctrl[dd] = P82_x;
	dd++;
    dctrl[dd] = P82_y;
	dd++;
	dctrl[dd] = P83;
	dd++;
	dctrl[dd] = P84;
	dd++;
	ictrl[ii] = P85;
	ii++;
    dctrl[dd] = P86_x;
	dd++;
    dctrl[dd] = P86_y;
	dd++;
	dctrl[dd] = P87;
	dd++;
	dctrl[dd] = P88;
	dd++;
	dctrl[dd] = P89_cm;
	dd++;
	dctrl[dd] = P89_cd;
	dd++;
	dctrl[dd] = P90;
	dd++;
	dctrl[dd] = P91;
	dd++;
    ictrl[ii] = P92;
	ii++;
    ictrl[ii] = P93;
	ii++;
	ictrl[ii] = P101;
	ii++;
	dctrl[dd] = P101_xm;
	dd++;
	dctrl[dd] = P101_ym;
	dd++;
	dctrl[dd] = P101_zs;
	dd++;
	dctrl[dd] = P101_ze;
	dd++;
	dctrl[dd] = P101_r1;
	dd++;
	dctrl[dd] = P101_r2;
	dd++;
	ictrl[ii] = P121;
	ii++;
	ictrl[ii] = P122;
	ii++;
	ictrl[ii] = P123;
	ii++;
	ictrl[ii] = P124;
	ii++;
	ictrl[ii] = P125;
	ii++;
	ictrl[ii] = P126;
	ii++;
	ictrl[ii] = P151;
	ii++;
	ictrl[ii] = P152;
	ii++;
	ictrl[ii] = P180;
	ii++;
	ictrl[ii] = P181;
	ii++;
	dctrl[dd] = P182;
	dd++;
    ictrl[ii] = P210;
	ii++;
	ictrl[ii] = P211;
	ii++;
	dctrl[dd] = P212;
	dd++;
    ictrl[ii] = P230;
	ii++;
	ictrl[ii] = P351;
	ii++;
	ictrl[ii] = P352;
	ii++;

    ictrl[ii] = S10;
	ii++;
    ictrl[ii] = S11;
	ii++;
    ictrl[ii] = S12;
	ii++;
    dctrl[dd] = S13;
	dd++;
    dctrl[dd] = S14;
	dd++;
    ictrl[ii] = S15;
	ii++;
    ictrl[ii] = S16;
	ii++;
    ictrl[ii] = S17;
	ii++;
    ictrl[ii] = S18;
	ii++;
	dctrl[dd] = S19;
	dd++;
    dctrl[dd] = S20;
	dd++;
	dctrl[dd] = S21;
	dd++;
    dctrl[dd] = S22;
	dd++;
    dctrl[dd] = S23;
	dd++;
    dctrl[dd] = S24;
	dd++;
    dctrl[dd] = S25;
	dd++;
	dctrl[dd] = S28;
	dd++;
	dctrl[dd] = S29;
	dd++;
    dctrl[dd] = S30;
	dd++;
	ictrl[ii] = S31;
	ii++;
    ictrl[ii] = S37;
	ii++;
	ictrl[ii] = S38;
	ii++;
	ictrl[ii] = S39;
	ii++;
    ictrl[ii] = S41;
	ii++;
	ictrl[ii] = S42;
	ii++;
    ictrl[ii] = S43;
	ii++;
    ictrl[ii] = S44;
	ii++;
	dctrl[dd] = S45;
	dd++;
	dctrl[dd] = S46;
	dd++;
	dctrl[dd] = S47;
	dd++;
	dctrl[dd] = S48;
	dd++;
    ictrl[ii] = S50;
	ii++;
    dctrl[dd] = S57;
	dd++;
    dctrl[dd] = S60;
	dd++;
    dctrl[dd] = S71;
	dd++;
    dctrl[dd] = S72;
	dd++;
    ictrl[ii] = S73;
	ii++;
    ictrl[ii] = S77;
	ii++;
    dctrl[dd] = S77_xs;
	dd++;
    dctrl[dd] = S77_xe;
	dd++;
    ictrl[ii] = S80;
	ii++;
    dctrl[dd] = S81;
	dd++;
    dctrl[dd] = S82;
	dd++;
    ictrl[ii] = S90;
	ii++;
    ictrl[ii] = S91;
	ii++;
	dctrl[dd] = S92;
	dd++;
    dctrl[dd] = S93;
	dd++;
	ictrl[ii] = S100;
	ii++;
	ictrl[ii] = S101;
	ii++;
    ictrl[ii] = S102;
	ii++;
    ictrl[ii] = S103;
	ii++;
    dctrl[dd] = S116;
	dd++;
    dctrl[dd] = S117;
	dd++;

    ictrl[ii] = T10;
    ii++;
    ictrl[ii] = T11;
    ii++;
    ictrl[ii] = T12;
    ii++;
	dctrl[dd] = T13;
    dd++;
    ictrl[ii] = T30;
    ii++;
    dctrl[dd] = T31;
    dd++;
	dctrl[dd] = T35;
    dd++;
	ictrl[ii] = T36;
    ii++;
	dctrl[dd] = T37;
    dd++;
    dctrl[dd] = T38;
    dd++;
    dctrl[dd] = T39;
    dd++;
    ictrl[ii] = T40;
    ii++;
	ictrl[ii] = T41;
    ii++;
	dctrl[dd] = T42;
    dd++;
	ictrl[ii] = T43;
    ii++;
    dctrl[dd] = T51;
    dd++;
	dctrl[dd] = T52;
    dd++;
	dctrl[dd] = T53;
    dd++;

    dctrl[dd] = W1;
    dd++;
    dctrl[dd] = W2;
    dd++;
    dctrl[dd] = W3;
    dd++;
    dctrl[dd] = W4;
    dd++;
    dctrl[dd] = W5;
    dd++;
	dctrl[dd] = W6;
    dd++;
	dctrl[dd] = W7;
    dd++;
    dctrl[dd] = W10;
    dd++;
    dctrl[dd] = W20;
    dd++;
    dctrl[dd] = W21;
    dd++;
    dctrl[dd] = W22;
    dd++;
	ictrl[ii] = W30;
    ii++;
	dctrl[dd] = W31;
    dd++;
    ictrl[ii] = W90;
    ii++;
	dctrl[dd] = W95;
    dd++;
    dctrl[dd] = W96;
    dd++;
    dctrl[dd] = W97;
    dd++;
	dctrl[dd] = W98;
    dd++;
    ictrl[ii] = W101;
    ii++;
	dctrl[dd] = W102_phi;
    dd++;
    dctrl[dd] = W102_c;
    dd++;
    dctrl[dd] = W103;
    dd++;
    dctrl[dd] = W104;
    dd++;
    ictrl[ii] = W110;
    ii++;
    ictrl[ii] = W111;
    ii++;
    dctrl[dd] = W112;
    dd++;
	
	ictrl[ii] = X10;
	ii++;
	ictrl[ii] = X11_u;
	ii++;
	ictrl[ii] = X11_v;
	ii++;
	ictrl[ii] = X11_w;
	ii++;
	ictrl[ii] = X11_p;
	ii++;
	ictrl[ii] = X11_q;
	ii++;
	ictrl[ii] = X11_r;
	ii++;
	ictrl[ii] = X12;
	ii++;
	ictrl[ii] = X13;
	ii++;  	
    ictrl[ii] = X18;
	ii++;
	ictrl[ii] = X19;
	ii++;
	ictrl[ii] = X21;
	ii++;
	dctrl[dd] = X21_d;
	dd++;
	ictrl[ii] = X22;
	ii++;
	dctrl[dd] = X22_m;
	dd++;
	ictrl[ii] = X23;
	ii++;
	dctrl[dd] = X23_x;
	dd++;
	dctrl[dd] = X23_y;
	dd++;
	dctrl[dd] = X23_z;
	dd++;
	ictrl[ii] = X24;
	ii++;
	dctrl[dd] = X24_Ix;
	dd++;
	dctrl[dd] = X24_Iy;
	dd++;
	dctrl[dd] = X24_Iz;
	dd++;
	dctrl[dd] = X25_Cp;
	dd++;
	dctrl[dd] = X25_Cq;
	dd++;
	dctrl[dd] = X25_Cr;
	dd++;
    dctrl[dd] = X26_Ku;
	dd++;
	dctrl[dd] = X26_Kv;
	dd++;
	dctrl[dd] = X26_Kw;
	dd++;
    ictrl[ii] = X27;
	ii++;
    dctrl[dd] = X27_x;
	dd++;
    dctrl[dd] = X27_y;
	dd++;
    dctrl[dd] = X27_z;
	dd++;
	ictrl[ii] = X31;
	ii++;
	ictrl[ii] = X32;
	ii++;
	ictrl[ii] = X33;
	ii++;
    ictrl[ii] = X34;
	ii++;
	ictrl[ii] = X38;
	ii++;
    ictrl[ii] = X40;
	ii++;
	dctrl[dd] = X41;
	dd++;
	ictrl[ii] = X100;
	ii++;
	dctrl[dd] = X100_x;
	dd++;
	dctrl[dd] = X100_y;
	dd++;
	dctrl[dd] = X100_z;
	dd++;
	ictrl[ii] = X101;
	ii++;
	dctrl[dd] = X101_phi;
	dd++;
	dctrl[dd] = X101_theta;
	dd++;
	dctrl[dd] = X101_psi;
	dd++;
	ictrl[ii] = X102;
	ii++;
	dctrl[dd] = X102_u;
	dd++;
	dctrl[dd] = X102_v;
	dd++;
	dctrl[dd] = X102_w;
	dd++;
	ictrl[ii] = X103;
	ii++;
	dctrl[dd] = X103_p;
	dd++;
	dctrl[dd] = X103_q;
	dd++;
	dctrl[dd] = X103_r;
	dd++;
	ictrl[ii] = X110;
	ii++;
	ictrl[ii] = X120;
	ii++;
	dctrl[dd] = X120_rad;
	dd++;
	dctrl[dd] = X120_xc;
	dd++;
	dctrl[dd] = X120_yc;
	dd++;
	dctrl[dd] = X120_zc;
	dd++;
	ictrl[ii] = X131;
	ii++;
	dctrl[dd] = X131_rad;
	dd++;
	dctrl[dd] = X131_h;
	dd++;
	dctrl[dd] = X131_xc;
	dd++;
	dctrl[dd] = X131_yc;
	dd++;
	dctrl[dd] = X131_zc;
	dd++;
	ictrl[ii] = X132;
	ii++;
	dctrl[dd] = X132_rad;
	dd++;
	dctrl[dd] = X132_h;
	dd++;
	dctrl[dd] = X132_xc;
	dd++;
	dctrl[dd] = X132_yc;
	dd++;
	dctrl[dd] = X132_zc;
	dd++;
	ictrl[ii] = X133;
	ii++;
	dctrl[dd] = X133_rad;
	dd++;
	dctrl[dd] = X133_h;
	dd++;
	dctrl[dd] = X133_xc;
	dd++;
	dctrl[dd] = X133_yc;
	dd++;
	dctrl[dd] = X133_zc;
	dd++;
	ictrl[ii] = X153;
	ii++;
	dctrl[dd] = X153_xs;
	dd++;
	dctrl[dd] = X153_xe;
	dd++;
	dctrl[dd] = X153_ys;
	dd++;
	dctrl[dd] = X153_ye;
	dd++;
	dctrl[dd] = X153_zs;
	dd++;
	dctrl[dd] = X153_ze;
	dd++;
    ictrl[ii] = X163;
	ii++;
    ictrl[ii] = X164;
	ii++;
    ictrl[ii] = X180;
	ii++;
    dctrl[dd] = X181;
	dd++;
    ictrl[ii] = X182;
	ii++;
    dctrl[dd] = X182_x;
	dd++;
    dctrl[dd] = X182_y;
	dd++;
    dctrl[dd] = X182_z;
	dd++;
    ictrl[ii] = X183;
	ii++;
    dctrl[dd] = X183_x;
	dd++;
    dctrl[dd] = X183_y;
	dd++;
    dctrl[dd] = X183_z;
	dd++;
    dctrl[dd] = X183_phi;
	dd++;
    dctrl[dd] = X183_theta;
	dd++;
    dctrl[dd] = X183_psi;
	dd++;
	ictrl[ii] = X210;
	ii++;
	dctrl[dd] = X210_u;
	dd++;
	dctrl[dd] = X210_v;
	dd++;
	dctrl[dd] = X210_w;
	dd++;
	ictrl[ii] = X211;
	ii++;
	dctrl[dd] = X211_p;
	dd++;
	dctrl[dd] = X211_q;
	dd++;
	dctrl[dd] = X211_r;
	dd++;
    ictrl[ii] = X221;
	ii++;
	dctrl[dd] = X221_xs;
	dd++;
    dctrl[dd] = X221_xe;
	dd++;
    dctrl[dd] = X221_ys;
	dd++;
    dctrl[dd] = X221_ye;
	dd++;
    dctrl[dd] = X221_zs;
	dd++;
    dctrl[dd] = X221_ze;
	dd++;
	ictrl[ii] = X310;
	ii++;		
    ictrl[ii] = X311;
	ii++;
    ictrl[ii] = X320;
	ii++;
	
	ictrl[ii] = Y50;
	ii++;
	ictrl[ii] = Y60;
	ii++;
    ictrl[ii] = Y71;
	ii++;
    ictrl[ii] = Y72;
	ii++;
    ictrl[ii] = Y73;
	ii++;
    ictrl[ii] = Y74;
	ii++;
    
    
    
	
// --------------------------

	
	for(n=0;n<B67;++n)
    {
    dctrl[dd] = B67_val[n];
    dd++;
	dctrl[dd] = B67_dist[n];
    dd++;
	dctrl[dd] = B67_b[n];
    dd++;
    dctrl[dd] = B67_x[n];
    dd++;
	dctrl[dd] = B67_y[n];
    dd++;
    }
	
	for(n=0;n<B68;++n)
    {
    dctrl[dd] = B68_val[n];
    dd++;
	dctrl[dd] = B68_dist[n];
    dd++;
	dctrl[dd] = B68_b[n];
    dd++;
    dctrl[dd] = B68_x[n];
    dd++;
	dctrl[dd] = B68_y[n];
    dd++;
    }
	
	for(n=0;n<B69;++n)
    {
    dctrl[dd] = B69_val[n];
    dd++;
	dctrl[dd] = B69_dist[n];
    dd++;
	dctrl[dd] = B69_b[n];
    dd++;
    dctrl[dd] = B69_x[n];
    dd++;
	dctrl[dd] = B69_y[n];
    dd++;
    }
	
	for(n=0;n<B70;++n)
    {
    dctrl[dd]   = B70_val[n];
    dd++;
	dctrl[dd] = B70_dist[n];
    dd++;
	dctrl[dd] = B70_b[n];
    dd++;
    dctrl[dd] = B70_x[n];
    dd++;
	dctrl[dd] = B70_y[n];
    dd++;
    }
	
	for(n=0;n<B71;++n)
    {
    dctrl[dd] = B71_val[n];
    dd++;
	dctrl[dd] = B71_dist[n];
    dd++;
	dctrl[dd] = B71_b[n];
    dd++;
    dctrl[dd] = B71_x[n];
    dd++;
	dctrl[dd] = B71_y[n];
    dd++;
    }
	
	for(n=0;n<B106;++n)
    {
    dctrl[dd] = B106_b[n];
    dd++;
    dctrl[dd] = B106_x[n];
    dd++;
	dctrl[dd] = B106_y[n];
    dd++;
    }
	
	for(n=0;n<B107;++n)
    {
    dctrl[dd] = B107_xs[n];
    dd++;
    dctrl[dd] = B107_xe[n];
    dd++;
	dctrl[dd] = B107_ys[n];
    dd++;
    dctrl[dd] = B107_ye[n];
    dd++;
    dctrl[dd] = B107_d[n];
    dd++;
    }
    
    for(n=0;n<B108;++n)
    {
    dctrl[dd] = B108_xs[n];
    dd++;
    dctrl[dd] = B108_xe[n];
    dd++;
	dctrl[dd] = B108_ys[n];
    dd++;
    dctrl[dd] = B108_ye[n];
    dd++;
    dctrl[dd] = B108_d[n];
    dd++;
    }
	
	for(n=0;n<B240;++n)
    {
    dctrl[dd] = B240_C[n];
    dd++;
    dctrl[dd] = B240_D[n];
    dd++;
	dctrl[dd] = B240_xs[n];
    dd++;
    dctrl[dd] = B240_xe[n];
    dd++;
    dctrl[dd] = B240_ys[n];
    dd++;
    dctrl[dd] = B240_ye[n];
    dd++;
    dctrl[dd] = B240_zs[n];
    dd++;
    dctrl[dd] = B240_ze[n];
    dd++;
    }
    
    for(n=0;n<B270;++n)
    {
	dctrl[dd] = B270_xs[n];
    dd++;
    dctrl[dd] = B270_xe[n];
    dd++;
    dctrl[dd] = B270_ys[n];
    dd++;
    dctrl[dd] = B270_ye[n];
    dd++;
    dctrl[dd] = B270_zs[n];
    dd++;
    dctrl[dd] = B270_ze[n];
    dd++;
    dctrl[dd] = B270_n[n];
    dd++;
    dctrl[dd] = B270_d50[n];
    dd++;
	dctrl[dd] = B270_alpha[n];
    dd++;
	dctrl[dd] = B270_beta[n];
    dd++;
    }
    
    for(n=0;n<B274;++n)
    {
	dctrl[dd] = B274_xc[n];
    dd++;
    dctrl[dd] = B274_yc[n];
    dd++;
    dctrl[dd] = B274_zs[n];
    dd++;
    dctrl[dd] = B274_ze[n];
    dd++;
    dctrl[dd] = B274_r[n];
    dd++;
    dctrl[dd] = B274_n[n];
    dd++;
    dctrl[dd] = B274_d50[n];
    dd++;
	dctrl[dd] = B274_alpha[n];
    dd++;
	dctrl[dd] = B274_beta[n];
    dd++;
    }
    
    for(n=0;n<B281;++n)
    {
	dctrl[dd] = B281_xs[n];
    dd++;
    dctrl[dd] = B281_xe[n];
    dd++;
    dctrl[dd] = B281_ys[n];
    dd++;
    dctrl[dd] = B281_ye[n];
    dd++;
    dctrl[dd] = B281_zs[n];
    dd++;
    dctrl[dd] = B281_ze[n];
    dd++;
    dctrl[dd] = B281_n[n];
    dd++;
    dctrl[dd] = B281_d50[n];
    dd++;
	dctrl[dd] = B281_alpha[n];
    dd++;
	dctrl[dd] = B281_beta[n];
    dd++;
    }
    
    for(n=0;n<B291;++n)
    {
	dctrl[dd] = B291_xs[n];
    dd++;
    dctrl[dd] = B291_xe[n];
    dd++;
    dctrl[dd] = B291_ys[n];
    dd++;
    dctrl[dd] = B291_ye[n];
    dd++;
    dctrl[dd] = B291_zs[n];
    dd++;
    dctrl[dd] = B291_ze[n];
    dd++;
    dctrl[dd] = B291_d[n];
    dd++;
    dctrl[dd] = B291_n[n];
    dd++;
    dctrl[dd] = B291_d50[n];
    dd++;
	dctrl[dd] = B291_alpha[n];
    dd++;
	dctrl[dd] = B291_beta[n];
    dd++;
    }
    
    for(n=0;n<B310;++n)
    {
	dctrl[dd] = B310_xs[n];
    dd++;
    dctrl[dd] = B310_xe[n];
    dd++;
    dctrl[dd] = B310_ys[n];
    dd++;
    dctrl[dd] = B310_ye[n];
    dd++;
    dctrl[dd] = B310_zs[n];
    dd++;
    dctrl[dd] = B310_ze[n];
    dd++;
    dctrl[dd] = B310_N[n];
    dd++;
    dctrl[dd] = B310_D[n];
    dd++;
	dctrl[dd] = B310_Cd[n];
    dd++;
    }
	
	for(n=0;n<C75;++n)
    {
    dctrl[dd] = C75_x[n];
    dd++;
    dctrl[dd] = C75_z[n];
    dd++;
    dctrl[dd] = C75_a[n];
    dd++;
    dctrl[dd] = C75_s[n];
    dd++;
    dctrl[dd] = C75_l[n];
    dd++;
	dctrl[dd] = C75_v[n];
    dd++;
    }
	
    for(n=0;n<F70;++n)
    {
    dctrl[dd] = F70_xs[n];
	dd++;
    dctrl[dd] = F70_xe[n];
	dd++;
    dctrl[dd] = F70_ys[n];
	dd++;
    dctrl[dd] = F70_ye[n];
	dd++;
    dctrl[dd] = F70_zs[n];
	dd++;
    dctrl[dd] = F70_ze[n];
	dd++;
    }
	
	for(n=0;n<F71;++n)
    {
    dctrl[dd] = F71_xs[n];
	dd++;
    dctrl[dd] = F71_xe[n];
	dd++;
    dctrl[dd] = F71_ys[n];
	dd++;
    dctrl[dd] = F71_ye[n];
	dd++;
    dctrl[dd] = F71_zs[n];
	dd++;
    dctrl[dd] = F71_ze[n];
	dd++;
    }
	
	for(n=0;n<F72;++n)
    {
    dctrl[dd] = F72_xs[n];
	dd++;
    dctrl[dd] = F72_xe[n];
	dd++;
    dctrl[dd] = F72_ys[n];
	dd++;
    dctrl[dd] = F72_ye[n];
	dd++;
    dctrl[dd] = F72_h[n];
	dd++;
    }
	
    for(n=0;n<F369;++n)
    {
    dctrl[dd] = F369_x[n];
	dd++;
    dctrl[dd] = F369_z[n];
	dd++;
    dctrl[dd] = F369_a[n];
	dd++;
    dctrl[dd] = F369_s[n];
	dd++;
    dctrl[dd] = F369_l[n];
	dd++;
	dctrl[dd] = F369_v[n];
	dd++;
    }
    
	for(n=0;n<F370;++n)
    {
    dctrl[dd] = F370_xs[n];
	dd++;
    dctrl[dd] = F370_xe[n];
	dd++;
    dctrl[dd] = F370_ys[n];
	dd++;
    dctrl[dd] = F370_ye[n];
	dd++;
    dctrl[dd] = F370_zs[n];
	dd++;
    dctrl[dd] = F370_ze[n];
	dd++;
    }
	
	for(n=0;n<F371;++n)
    {
    dctrl[dd] = F371_xs[n];
	dd++;
    dctrl[dd] = F371_xe[n];
	dd++;
    dctrl[dd] = F371_ys[n];
	dd++;
    dctrl[dd] = F371_ye[n];
	dd++;
    dctrl[dd] = F371_zs[n];
	dd++;
    dctrl[dd] = F371_ze[n];
	dd++;
    }
    
    for(n=0;n<F374;++n)
    {
    dctrl[dd] = F374_xc[n];
	dd++;
    dctrl[dd] = F374_zc[n];
	dd++;
    dctrl[dd] = F374_r[n];
	dd++;
    }
    
    for(n=0;n<F375;++n)
    {
    dctrl[dd] = F375_xc[n];
	dd++;
    dctrl[dd] = F375_zc[n];
	dd++;
    dctrl[dd] = F375_r[n];
	dd++;
    }
    
    for(n=0;n<F378;++n)
    {
    dctrl[dd] = F378_xc[n];
	dd++;
    dctrl[dd] = F378_yc[n];
	dd++;
    dctrl[dd] = F378_zc[n];
	dd++;
    dctrl[dd] = F378_r[n];
	dd++;
    }
    
    for(n=0;n<F379;++n)
    {
    dctrl[dd] = F379_xc[n];
	dd++;
    dctrl[dd] = F379_yc[n];
	dd++;
    dctrl[dd] = F379_zc[n];
	dd++;
    dctrl[dd] = F379_r[n];
	dd++;
    }
	
	for(n=0;n<F390;++n)
    {
    dctrl[dd] = F390_xs[n];
	dd++;
    dctrl[dd] = F390_xe[n];
	dd++;
    dctrl[dd] = F390_ys[n];
	dd++;
    dctrl[dd] = F390_ye[n];
	dd++;
    dctrl[dd] = F390_zs[n];
	dd++;
    dctrl[dd] = F390_ze[n];
	dd++;
    }
	
	for(n=0;n<F391;++n)
    {
    dctrl[dd] = F391_xs[n];
	dd++;
    dctrl[dd] = F391_xe[n];
	dd++;
    dctrl[dd] = F391_ys[n];
	dd++;
    dctrl[dd] = F391_ye[n];
	dd++;
    dctrl[dd] = F391_zs[n];
	dd++;
    dctrl[dd] = F391_ze[n];
	dd++;
    }
    
    for(n=0;n<F394;++n)
    {
    dctrl[dd] = F394_xc[n];
	dd++;
    dctrl[dd] = F394_zc[n];
	dd++;
    dctrl[dd] = F394_r[n];
	dd++;
    }
    
    for(n=0;n<F395;++n)
    {
    dctrl[dd]   = F395_xc[n];
	dd++;
    dctrl[dd] = F395_zc[n];
	dd++;
    dctrl[dd] = F395_r[n];
	dd++;
    }
    
    for(n=0;n<F398;++n)
    {
    dctrl[dd] = F398_xc[n];
	dd++;
    dctrl[dd] = F398_yc[n];
	dd++;
    dctrl[dd] = F398_zc[n];
	dd++;
    dctrl[dd] = F398_r[n];
	dd++;
    }
    
    for(n=0;n<F399;++n)
    {
    dctrl[dd] = F399_xc[n];
	dd++;
    dctrl[dd] = F399_yc[n];
	dd++;
    dctrl[dd] = F399_zc[n];
	dd++;
    dctrl[dd] = F399_r[n];
	dd++;
    }

    for(n=0;n<G60;++n)
    {
    dctrl[dd] = G60_xs[n];
	dd++;
    dctrl[dd] = G60_xe[n];
	dd++;
    dctrl[dd] = G60_ys[n];
	dd++;
    dctrl[dd] = G60_ye[n];
	dd++;
    dctrl[dd] = G60_zs[n];
	dd++;
    dctrl[dd] = G60_ze[n];
	dd++;
    }

    for(n=0;n<G61;++n)
    {
    dctrl[dd] = G61_xs[n];
	dd++;
    dctrl[dd] = G61_xe[n];
	dd++;
    dctrl[dd] = G61_ys[n];
	dd++;
    dctrl[dd] = G61_ye[n];
	dd++;
    dctrl[dd] = G61_zs[n];
	dd++;
    dctrl[dd] = G61_ze[n];
	dd++;
    }
	
	for(n=0;n<G81;++n)
    {
    dctrl[dd] = G81_xs[n];
    dd++;
    dctrl[dd] = G81_xe[n];
    dd++;
    dctrl[dd] = G81_ys[n];
    dd++;
    dctrl[dd] = G81_ye[n];
    dd++;
    dctrl[dd] = G81_zs[n];
    dd++;
    dctrl[dd] = G81_ze[n];
    dd++;
    }
	
	for(n=0;n<G95;++n)
    {
    dctrl[dd] = G95_xs[n];
    dd++;
    dctrl[dd] = G95_xe[n];
    dd++;
    dctrl[dd] = G95_ys[n];
    dd++;
    dctrl[dd] = G95_ye[n];
    dd++;
    dctrl[dd] = G95_zs[n];
    dd++;
    dctrl[dd] = G95_ze[n];
    dd++;
    }
	
	for(n=0;n<P35;++n)
    {
    dctrl[dd]   = P35_ts[n];
	dd++;
	dctrl[dd]  = P35_te[n];
	dd++;
    dctrl[dd] = P35_dt[n];
	dd++;
    }

    for(n=0;n<P50;++n)
    {
    dctrl[dd] = P50_x[n];
	dd++;
    dctrl[dd] = P50_y[n];
	dd++;
    }
	
	for(n=0;n<P51;++n)
    {
    dctrl[dd] = P51_x[n];
	dd++;
    dctrl[dd] = P51_y[n];
	dd++;
    }

    for(n=0;n<P52;++n)
	{
    dctrl[dd] = P52_y[n];
	dd++;
	}
	
	for(n=0;n<P56;++n)
	{
    dctrl[dd] = P56_x[n];
	dd++;
	}
	
	for(n=0;n<P61;++n)
    {
    dctrl[dd] = P61_x[n];
	dd++;
    dctrl[dd] = P61_y[n];
	dd++;
	dctrl[dd] = P61_z[n];
	dd++;
    }
	
	for(n=0;n<P62;++n)
    {
    dctrl[dd] = P62_xs[n];
	dd++;
    dctrl[dd] = P62_xe[n];
	dd++;
	dctrl[dd] = P62_ys[n];
	dd++;
	dctrl[dd] = P62_ye[n];
	dd++;
    dctrl[dd] = P62_zs[n];
	dd++;
	dctrl[dd] = P62_ze[n];
	dd++;
    }
	
	for(n=0;n<P67;++n)
	{
    dctrl[dd] = P67_x[n];
	dd++;
	}
	
	for(n=0;n<P81;++n)
    {
    dctrl[dd] = P81_xs[n];
	dd++;
    dctrl[dd] = P81_xe[n];
	dd++;
	dctrl[dd] = P81_ys[n];
	dd++;
	dctrl[dd] = P81_ye[n];
	dd++;
    dctrl[dd] = P81_zs[n];
	dd++;
	dctrl[dd] = P81_ze[n];
	dd++;
    }
	
	for(n=0;n<P85;++n)
    {
    dctrl[dd] = P85_xs[n];
	dd++;
    dctrl[dd] = P85_xe[n];
	dd++;
	dctrl[dd] = P85_ys[n];
	dd++;
	dctrl[dd] = P85_ye[n];
	dd++;
    dctrl[dd] = P85_zs[n];
	dd++;
	dctrl[dd] = P85_ze[n];
	dd++;
    }
	
	for(n=0;n<P121;++n)
    {
    dctrl[dd] = P121_x[n];
	dd++;
    dctrl[dd] = P121_y[n];
	dd++;
    }
	
	for(n=0;n<P123;++n)
	{
    dctrl[dd] = P123_y[n];
	dd++;
	}
	
	for(n=0;n<P124;++n)
	{
    dctrl[dd] = P124_x[n];
	dd++;
	}
	
	for(n=0;n<P125;++n)
    {
    dctrl[dd] = P125_x[n];
	dd++;
    dctrl[dd] = P125_y[n];
	dd++;
    }
    
    for(n=0;n<P230;++n)
    {
    dctrl[dd] = P230_x[n];
	dd++;
    }
	
	for(n=0;n<P351;++n)
    {
    dctrl[dd] = P351_x[n];
	dd++;
    dctrl[dd] = P351_y[n];
	dd++;
    }
	
	for(n=0;n<P352;++n)
    {
    dctrl[dd] = P352_x[n];
	dd++;
    dctrl[dd] = P352_y[n];
	dd++;
    }
    
    for(n=0;n<S73;++n)
    {
    dctrl[dd] = S73_val[n];
    dd++;
	dctrl[dd] = S73_dist[n];
    dd++;
	dctrl[dd] = S73_b[n];
    dd++;
    dctrl[dd] = S73_x[n];
    dd++;
	dctrl[dd] = S73_y[n];
    dd++;
    }
	
	for(n=0;n<X110;++n)
    {
    dctrl[dd] = X110_xs[n];
    dd++;
    dctrl[dd] = X110_xe[n];
    dd++;
    dctrl[dd] = X110_ys[n];
    dd++;
    dctrl[dd] = X110_ye[n];
    dd++;
    dctrl[dd] = X110_zs[n];
    dd++;
    dctrl[dd] = X110_ze[n];
    dd++;
    }
    
    for(n=0;n<X163;++n)
    {
    dctrl[dd] = X163_x1[n];
    dd++;
    dctrl[dd] = X163_y1[n];
    dd++;
    dctrl[dd] = X163_z1[n];
    dd++;
    dctrl[dd] = X163_x2[n];
    dd++;
    dctrl[dd] = X163_y2[n];
    dd++;
    dctrl[dd] = X163_z2[n];
    dd++;
    dctrl[dd] = X163_x3[n];
    dd++;
    dctrl[dd] = X163_y3[n];
    dd++;
    dctrl[dd] = X163_z3[n];
    dd++;
    dctrl[dd] = X163_x4[n];
    dd++;
    dctrl[dd] = X163_y4[n];
    dd++;
    dctrl[dd] = X163_z4[n];
    dd++;
    dctrl[dd] = X163_x5[n];
    dd++;
    dctrl[dd] = X163_y5[n];
    dd++;
    dctrl[dd] = X163_z5[n];
    dd++;
    dctrl[dd] = X163_x6[n];
    dd++;
    dctrl[dd] = X163_y6[n];
    dd++;
    dctrl[dd] = X163_z6[n];
    dd++;
    }
    
    for(n=0;n<X164;++n)
    {
    dctrl[dd] = X164_x1[n];
    dd++;
    dctrl[dd] = X164_y1[n];
    dd++;
    dctrl[dd] = X164_z1[n];
    dd++;
    dctrl[dd] = X164_x2[n];
    dd++;
    dctrl[dd] = X164_y2[n];
    dd++;
    dctrl[dd] = X164_z2[n];
    dd++;
    dctrl[dd] = X164_x3[n];
    dd++;
    dctrl[dd] = X164_y3[n];
    dd++;
    dctrl[dd] = X164_z3[n];
    dd++;
    dctrl[dd] = X164_x4[n];
    dd++;
    dctrl[dd] = X164_y4[n];
    dd++;
    dctrl[dd] = X164_z4[n];
    dd++;
    dctrl[dd] = X164_x5[n];
    dd++;
    dctrl[dd] = X164_y5[n];
    dd++;
    dctrl[dd] = X164_z5[n];
    dd++;
    dctrl[dd] = X164_x6[n];
    dd++;
    dctrl[dd] = X164_y6[n];
    dd++;
    dctrl[dd] = X164_z6[n];
    dd++;
    dctrl[dd] = X164_x7[n];
    dd++;
    dctrl[dd] = X164_y7[n];
    dd++;
    dctrl[dd] = X164_z7[n];
    dd++;
    dctrl[dd] = X164_x8[n];
    dd++;
    dctrl[dd] = X164_y8[n];
    dd++;
    dctrl[dd] = X164_z8[n];
    dd++;
    }
    
    for(n=0;n<X311;++n)
    {
    dctrl[dd] = X311_xs[n];
    dd++;
    dctrl[dd] = X311_xe[n];
    dd++;
    dctrl[dd] = X311_ys[n];
    dd++;
    dctrl[dd] = X311_ye[n];
    dd++;
    dctrl[dd] = X311_zs[n];
    dd++;
    dctrl[dd] = X311_ze[n];
    dd++;
	dctrl[dd] = X311_w[n];
    dd++;
    dctrl[dd] = X311_rho_c[n];
    dd++;
    dctrl[dd] = X311_EA[n];
    dd++;
    dctrl[dd] = X311_d[n];
    dd++;
    dctrl[dd] = X311_l[n];
    dd++;
    dctrl[dd] = X311_H[n];
    dd++; 
    dctrl[dd] = X311_P[n];
    dd++;  
    dctrl[dd] = X311_facT[n];
    dd++;   
    }
}



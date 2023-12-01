/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"lexer.h"
#include <fstream>
#include <ctype.h>

void lexer::read_control()
{
    std::string line;
	char c;
	int numint;
    int count=0;

	ifstream control("ctrl.txt", ios_base::in);
	if(!control)
	{
		cout<<endl<<("no 'ctrl.txt' file found")<<endl<<endl;
		exit(0);
	}

	if(mpirank==0)
	cout<<"read ctrl"<<endl;


	while(!control.eof())
	{
		control>>c;

	if (c == '/')
	{
	control.ignore(1000, '\n');
	}
	else
	{
		switch(c)
		{
        case 'A': control>>numint;
				switch(numint)
				{
				case 10: control>>A10;
						 clear(c,numint);
						 break;
				case 210: control>>A210;
						 clear(c,numint);
						 break;
               case 209: control>>A209;
						 clear(c,numint);
						 break;
				case 211: control>>A211;
						 clear(c,numint);
						 break;
				case 212: control>>A212;
						 clear(c,numint);
						 break;
               case 214: control>>A214;
						 clear(c,numint);
						 break;
               case 215: control>>A215;
						 clear(c,numint);
						 break;
               case 216: control>>A216;
						 clear(c,numint);
						 break;
               case 217: control>>A217;
						 clear(c,numint);
						 break;
               case 218: control>>A218;
						 clear(c,numint);
						 break;
               case 219: control>>A219;
						 clear(c,numint);
						 break;
				case 220: control>>A220;
						 clear(c,numint);
						 break;
               case 221: control>>A221;
						 clear(c,numint);
						 break;
               case 223: control>>A223;
						 clear(c,numint);
						 break;
				case 230: control>>A230;
						 clear(c,numint);
						 break;
               case 240: control>>A240;
						 clear(c,numint);
						 break;
				case 241: control>>A241;
						 clear(c,numint);
						 break;
				case 242: control>>A242;
						 clear(c,numint);
						 break;
               case 243: control>>A243;
						 clear(c,numint);
						 break;
               case 244: control>>A244_val;
                        A244=1;
						 clear(c,numint);
						 break;
               case 245: control>>A245_val;
                        A245=1;
						 clear(c,numint);
						 break;
               case 246: control>>A246;
						 clear(c,numint);
						 break;
               case 247: control>>A247;
                        clear(c,numint);
						 break;
               case 248: control>>A248;
						 clear(c,numint);
						 break;
               case 249: control>>A249;
						 clear(c,numint);
						 break;
               case 250: control>>A250;
						 clear(c,numint);
						 break;
               case 251: control>>A251_val;
                        A251=1;
						 clear(c,numint);
						 break;
               case 260: control>>A260;
						 clear(c,numint);
						 break;
               case 261: control>>A261;
						 clear(c,numint);
						 break;
               case 262: control>>A262;
						 clear(c,numint);
						 break;
               case 310: control>>A310;
						 clear(c,numint);
						 break;
               case 311: control>>A311;
						 clear(c,numint);
						 break;
               case 312: control>>A312;
						 clear(c,numint);
						 break;
               case 313: control>>A313;
						 clear(c,numint);
						 break;
               case 320: control>>A320;
						 clear(c,numint);
						 break;
               case 321: control>>A321;
						 clear(c,numint);
						 break;
               case 322: control>>A322;
						 clear(c,numint);
						 break;
               case 323: control>>A323;
						 clear(c,numint);
						 break;
               case 329: control>>A329;
						 clear(c,numint);
						 break;
               case 340: control>>A340;
						 clear(c,numint);
						 break;
               case 341: control>>A341;
						 clear(c,numint);
						 break;
               case 342: control>>A342;
						 clear(c,numint);
						 break;
               case 343: control>>A343;
						 clear(c,numint);
						 break;
               case 344: control>>A344_val;
                        A344=1;
						 clear(c,numint);
						 break;
               case 345: control>>A345_val;
                        A345=1;
						 clear(c,numint);
						 break;
               case 346: control>>A346;
						 clear(c,numint);
						 break;
               case 347: control>>A347;
						 clear(c,numint);
						 break;
               case 348: control>>A348;
						 clear(c,numint);
						 break;
               case 350: control>>A350;
						 clear(c,numint);
						 break;
               case 351: control>>A351;
						 clear(c,numint);
						 break;
               case 352: control>>A352;
						 clear(c,numint);
						 break;
               case 353: control>>A353;
						 clear(c,numint);
						 break;
               case 354: control>>A354;
						 clear(c,numint);
						 break;
               case 355: control>>A355;
						 clear(c,numint);
						 break;
               case 356: control>>A356;
						 clear(c,numint);
						 break;
               case 357: control>>A357;
						 clear(c,numint);
						 break;
               case 361: control>>A361;
						 clear(c,numint);
						 break;
               case 362: control>>A362;
						 clear(c,numint);
						 break;
               case 363: control>>A363;
						 clear(c,numint);
						 break;
               case 365: control>>A365;
						 clear(c,numint);
						 break;
               case 368: control>>A368;
						 clear(c,numint);
						 break;
               case 369: control>>A369;
						 clear(c,numint);
						 break;
               case 410: control>>A410;
                        clear(c,numint);
                        break;
               case 440: control>>A440;
                        clear(c,numint);
                        break;
                        
               case 501: control>>A501;
                        clear(c,numint);
                        break;
               case 510: control>>A510;
                        clear(c,numint);
                        break;
               case 511: control>>A511;
                        clear(c,numint);
                        break;
               case 512: control>>A512;
                        clear(c,numint);
                        break;
               case 514: control>>A514;
                        clear(c,numint);
                        break;
               case 515: control>>A515;
                        clear(c,numint);
                        break;
               case 516: control>>A516;
                        clear(c,numint);
                        break;
               case 517: control>>A517;
                        clear(c,numint);
                        break;
               case 518: control>>A518;
                        clear(c,numint);
                        break;
               case 520: control>>A520;
                        clear(c,numint);
                        break;
               case 521: control>>A521;
                        clear(c,numint);
                        break;
               case 523: control>>A523;
                        clear(c,numint);
                        break;
               case 531: control>>A531;
                        clear(c,numint);
                        break;
               case 540: control>>A540;
                        clear(c,numint);
                        break;
               case 541: control>>A541;
                        clear(c,numint);
                        break;
               case 542: control>>A542;
                        clear(c,numint);
                        break;
               case 543: control>>A543;
                        clear(c,numint);
                        break;
               case 544: control>>A544;
                        clear(c,numint);
                        break;
               case 550: control>>A550;
                        clear(c,numint);
                        break;
               case 551: control>>A551;
                        clear(c,numint);
                        break;
               case 552: control>>A552;
                        clear(c,numint);
                        break;
               case 553: control>>A553;
                        clear(c,numint);
                        break;
				}
				break;
		case 'B': control>>numint;
				switch(numint)
				{
				case 10: control>>B10;
						 clear(c,numint);
						 break;
				case 20: control>>B20;
						 clear(c,numint);
						 break;
                 case 23: control>>B23;
						 clear(c,numint);
						 break;
				case 29: control>>B29;
						 clear(c,numint);
						 break;
				case 30: control>>B30;
						 clear(c,numint);
						 break;
                 case 31: control>>B31;
						 clear(c,numint);
						 break;
                 case 32: control>>B32_x>>B32_y>>B32_z;
                           B32=1;
						 clear(c,numint);
						 break;
                 case 33: control>>B33;
						 clear(c,numint);
						 break;
				case 50: control>>B50;
						 clear(c,numint);
						 break;
				case 51: control>>B51;
						 clear(c,numint);
						 break;
				case 52: control>>B52;
						 clear(c,numint);
						 break;
				case 53: control>>B53;
						 clear(c,numint);
						 break;
				case 54: control>>B54;
						 clear(c,numint);
						 break;
				case 55: control>>B55;
						 clear(c,numint);
						 break;
				case 56: control>>B56;
						 clear(c,numint);
						 break;
                case 60: control>>B60;
						 clear(c,numint);
						 break;
                case 61: control>>B61;
						 clear(c,numint);
						 break;
			   case 70: ++B70;
						 clear(c,numint);
						 break;
			   case 71: ++B71;
						 clear(c,numint);
						 break;
			   case 75: control>>B75;
						 clear(c,numint);
						 break;
              case 76: control>>B76;
						 clear(c,numint);
						 break;
			   case 77: control>>B77;
						 clear(c,numint);
						 break;
			   case 81: control>>B81_1>>B81_3>>B81_2;
						 B81=1;
						 clear(c,numint);
						 break;
              case 82: control>>B82;
						 clear(c,numint);
						 break;
			   case 83: control>>B83;
						 clear(c,numint);
						 break;
			   case 84: control>>B84;
						 clear(c,numint);
						 break;
              case 85: control>>B85;
						 clear(c,numint);
						 break;
			   case 86: control>>B86;
						 clear(c,numint);
						 break;
			   case 87: control>>B87_1>>B87_2;
                         B87=1;
						 clear(c,numint);
						 break;
			   case 88: control>>B88;
						 clear(c,numint);
						 break;
				case 89: control>>B89;
						 clear(c,numint);
						 break;
               case 90: control>>B90;
						 clear(c,numint);
						 break;
               case 91: {
                            control>>B91_1>>B91_2;
                            
                            int pos = control.tellg();
                            string test;
                            control>>test;
                            if (!isdigit(test[0])) 
                            {
                                control.seekg(pos);
                            }
                            else
                            {
                                cout<<endl;
                                cout<<"!!! wrong input error for B 91 !!!"<<endl<<endl;
                                cout<<"!!! please check the REEF3D User Guide !!!"<<endl<<endl<<endl<<endl;
                                exit(0);
                            }
                            B91=1;
						    clear(c,numint);
						    break;
                        }
               case 92: control>>B92;
						 clear(c,numint);
						 break;
               case 93: {
                            control>>B93_1>>B93_2;
                            
                            int pos = control.tellg();
                            string test;
                            control>>test;
                            if (!isdigit(test[0])) 
                            {
                                control.seekg(pos);
                            }
                            else
                            {
                                cout<<endl;
                                cout<<"!!! wrong input error for B 93 !!!"<<endl<<endl;
                                cout<<"!!! please check the REEF3D User Guide !!!"<<endl<<endl<<endl<<endl;
                                exit(0);
                            }
                            B93=1;
						    clear(c,numint);
						    break;
                        }
              case 94: control>>B94_wdt;
                        B94=1;
						 clear(c,numint);
						 break;
               case 96: {
                            control>>B96_1>>B96_2;
                            
                            int pos = control.tellg();
                            string test;
                            control>>test;
                            if (!isdigit(test[0])) 
                            {
                                control.seekg(pos);
                            }
                            else
                            {
                                cout<<endl;
                                cout<<"!!! wrong input error for B 96 !!!"<<endl<<endl;
                                cout<<"!!! please check the REEF3D User Guide !!!"<<endl<<endl<<endl<<endl;
                                exit(0);
                            }
						    clear(c,numint);
						    break;
                        }
               case 98: control>>B98;
						 clear(c,numint);
						 break;
               case 99: control>>B99;
						 clear(c,numint);
						 break;
               case 101: control>>B101;
						 clear(c,numint);
						 break;
               case 102: control>>B102;
						 clear(c,numint);
						 break;
			   case 105: control>>B105_1>>B105_2>>B105_3;
                        B105=1;
						 clear(c,numint);
						 break;
			   case 106: ++B106;
						 clear(c,numint);
						 break;
			   case 107: ++B107;
						 clear(c,numint);
						 break;
              case 108: ++B108;
						 clear(c,numint);
						 break;
              case 110: control>>B110_zs>>B110_ze;
                        B110=1;
						 clear(c,numint);
						 break;
			   case 111: control>>B111_zs>>B111_ze;
						 clear(c,numint);
						 break;
               case 112: control>>B112_zs>>B112_z2>>B112_ze;
						 clear(c,numint);
						 break;
               case 115: control>>B115;
						 clear(c,numint);
						 break;
               case 116: control>>B116;
						 clear(c,numint);
						 break;
               case 117: control>>B117;
						 clear(c,numint);
						 break;
			   case 120: control>>B120;
						 clear(c,numint);
						 break;
               case 122: control>>B122;
						 clear(c,numint);
						 break;
               case 123: control>>B123;
						 clear(c,numint);
						 break;
               case 125: control>>B125_y;
                           B125=1;
						 clear(c,numint);
						 break;
               case 127: control>>B127;
						 clear(c,numint);
						 break;
               case 130: control>>B130;
						 clear(c,numint);
						 break;
               case 131: control>>B131;
						 clear(c,numint);
						 break;
               case 132: control>>B132_s>>B132_e;
						 clear(c,numint);
						 break;
               case 133: control>>B133;
						 clear(c,numint);
						 break;
               case 134: control>>B134;
						 clear(c,numint);
						 break;
               case 135: control>>B135;
						 clear(c,numint);
						 break;
               case 136: control>>B136;
						 clear(c,numint);
						 break;
               case 138: control>>B138_1>>B138_2;
                           B138=1;
						 clear(c,numint);
						 break;
               case 139: control>>B139;
						 clear(c,numint);
						 break;
               case 140: control>>B140_1>>B140_2>>B140_3;
						 clear(c,numint);
						 break;
               case 160: control>>B160;
						 clear(c,numint);
						 break;
               case 170: control>>B170;
						 clear(c,numint);
						 break;
               case 180: control>>B180;
						 clear(c,numint);
						 break;
               case 181: control>>B181_1>>B181_2>>B181_3;
						 clear(c,numint);
                           B181=1;
						 break;
               case 182: control>>B182_1>>B182_2>>B182_3;
						 clear(c,numint);
                           B182=1;
						 break;
               case 183: control>>B183_1>>B183_2>>B183_3;
						 clear(c,numint);
                           B183=1;
						 break;
			   case 191: control>>B191_1>>B191_2>>B191_3>>B191_4;
						 clear(c,numint);
						 B191=1;
						 break;
			   case 192: control>>B192_1>>B192_2>>B192_3>>B192_4;
						 clear(c,numint);
						 B192=1;
						 break;
				case 194: control>>B194_s>>B194_e;
						 clear(c,numint);
						 break;
			   case 240: ++B240;
						 clear(c,numint);
						 break;
			   case 241: control>>B241;
						 clear(c,numint);
						 break;
			   case 242: control>>B242;
						 clear(c,numint);
						 break;
			   case 243: control>>B243;
						 clear(c,numint);
						 break;
               case 260: control>>B260;
						 clear(c,numint);
						 break;
               case 264: control>>B264;
						 clear(c,numint);
						 break;
               case 267: control>>B267;
						 clear(c,numint);
						 break;
               case 270: ++B270;
						 B269=1;
						 clear(c,numint);
						 break;
               case 274: ++B274;
						 B269=1;
						 clear(c,numint);
						 break;
               case 281: ++B281;
						 B269=1;
						 clear(c,numint);
						 break;
               case 282: ++B282;
						 B269=1;
						 clear(c,numint);
						 break;
               case 291: ++B291;
						 B269=1;
						 clear(c,numint);
						 break;
               case 295: control>>B295;
						 clear(c,numint);
						 break;
               case 308: control>>B308;
						 B269=2;
						 clear(c,numint);
						 break;
               case 309: control>>B309;
						 clear(c,numint);
						 break;
               case 310: ++B310;
						 clear(c,numint);
						 break;
               case 321: ++B321;
						 clear(c,numint);
						 break;
               case 322: ++B322;
						 clear(c,numint);
						 break;
               case 411: ++B411;
						 clear(c,numint);
						 break;
               case 412: ++B412;
						 clear(c,numint);
						 break;
               case 413: ++B413;
						 clear(c,numint);
						 break;
               case 414: ++B414;
						 clear(c,numint);
						 break;
               case 415: ++B415;
						 clear(c,numint);
						 break;
               case 416: ++B416;
						 clear(c,numint);
						 break;
               case 417: ++B417;
						 clear(c,numint);
						 break;
               case 418: ++B418;
						 clear(c,numint);
						 break;
               case 421: ++B421;
						 clear(c,numint);
						 break;
               case 422: ++B422;
						 clear(c,numint);
						 break;
               case 440: ++B440;
						 clear(c,numint);
						 break;
               case 441: ++B441;
						 clear(c,numint);
						 break;
               case 442: ++B442;
						 clear(c,numint);
						 break;
				}
				break;

		case 'C': control>>numint;
				switch(numint)
				{
				case 1: control>>C1;
						 clear(c,numint);
						 break;
				case 2: control>>C2;
						 clear(c,numint);
						 break;
				case 3: control>>C3;
						 clear(c,numint);
						 break;
				case 4: control>>C4;
						 clear(c,numint);
						 break;
				case 5: control>>C5;
						 clear(c,numint);
						 break;
               case 9: control>>C9;
						 clear(c,numint);
						 break;
				case 10: control>>C10;
						 clear(c,numint);
						 break;
				case 15: control>>C15;
						 clear(c,numint);
						 break;
				case 20: control>>C20;
						 clear(c,numint);
						 break;
                case 50: control>>C50_1>>C50_2;
						 clear(c,numint);
						 break;
                case 51: control>>C51;
						 clear(c,numint);
						 break;
                case 52: control>>C52;
						 clear(c,numint);
						 break;
                case 53: control>>C53;
						 clear(c,numint);
						 break;
                case 54: control>>C54;
						 clear(c,numint);
						 break;
                case 55: control>>C55;
						 clear(c,numint);
						 break;
                case 56: control>>C56;
						 clear(c,numint);
						 break;
                case 57: control>>C57_1>>C57_2>>C57_3>>C57_4;
						 clear(c,numint);
						 break;
                case 58: control>>C58_1>>C58_2>>C58_3>>C58_4;
						 clear(c,numint);
						 break;
				case 75: ++C75;
						 clear(c,numint);
						 break;

				}
				break;



		case 'D':control>>numint;
				switch(numint)
				{
				case 10: control>>D10;
						 clear(c,numint);
						 break;
				case 11: control>>D11;
						 clear(c,numint);
						 break;
				case 20: control>>D20;
						 clear(c,numint);
						 break;
				case 21: control>>D21;
						 clear(c,numint);
						 break;
				case 30: control>>D30;
						 clear(c,numint);
						 break;
                case 31: control>>D31;
						 clear(c,numint);
						 break;
                case 32: control>>D32;
						 clear(c,numint);
						 break;
                case 33: control>>D33;
						 clear(c,numint);
						 break;
                case 37: control>>D37;
						 clear(c,numint);
						 break;       
				}
				break;

		case 'F': control>>numint;
				switch(numint)
				{
				case 10: control>>F10;
						 clear(c,numint);
						 break;
				case 30: control>>F30;
						 clear(c,numint);
						 break;
                case 31: control>>F31;
						 clear(c,numint);
						 break;
                case 32: control>>F32;
						 clear(c,numint);
						 break;
                case 33: control>>F33;
						 clear(c,numint);
						 break;
                case 34: control>>F34;
						 clear(c,numint);
						 break;
                case 35: control>>F35;
						 clear(c,numint);
						 break;
				case 36: control>>F36;
						 clear(c,numint);
						 break;
				case 39: control>>F39;
						 clear(c,numint);
						 break;
				case 40: control>>F40;
						 clear(c,numint);
						 break;
				case 42: control>>F42;
						 clear(c,numint);
						 break;
				case 43: control>>F43;
						 clear(c,numint);
						 break;
				case 44: control>>F44;
						 clear(c,numint);
						 break;
                case 45: control>>F45;
						 clear(c,numint);
						 break;
                case 46: control>>F46;
						 clear(c,numint);
						 break;
                case 47: control>>F47;
						 clear(c,numint);
						 break;
				case 49: control>>F49;
						 clear(c,numint);
						 break;
				case 50: control>>F50;
						 clear(c,numint);
						 break;
				case 51: control>>F51;
                        F50_flag=1;
						 clear(c,numint);
						 break;
                case 52: control>>F52;
                        F50_flag=1;
						 clear(c,numint);
						 break;
                case 53: control>>F53;
                        F50_flag=1;
						 clear(c,numint);
						 break;
                case 54: control>>F54;
                        F50_flag=1;
						 clear(c,numint);
						 break;
                case 55: control>>F55;
                        F50_flag=1;
						 clear(c,numint);
						 break;
                case 56: control>>F56;
                        F50_flag=1;
						 clear(c,numint);
						 break;
                case 57: control>>F57_1>>F57_2>>F57_3>>F57_4;
						 clear(c,numint);
						 break;
                case 58: control>>F58_1>>F58_2>>F58_3>>F58_4;
						 clear(c,numint);
						 break;
                case 59: control>>F59_xm>>F59_ym>>F59_zs>>F59_ze>>F59_r;
						 clear(c,numint);
						 break;
                case 60: control>>F60;
						 clear(c,numint);
						 break;
                case 61: control>>F61;
						 clear(c,numint);
						 break;
                case 62: control>>F62;
						 clear(c,numint);
						 break;
                case 63: control>>F63;
						 clear(c,numint);
						 break;
                case 64: control>>F64_xs>>F64_ys>>F64_zs>>F64_alpha;
						 F64=1;
						 clear(c,numint);
						 break;
                case 70: ++F70;
						 clear(c,numint);
						 break;
				case 71: ++F71;
						 clear(c,numint);
						 break;
				case 72: ++F72;
						 clear(c,numint);
						 break;
                case 80: control>>F80;
						 clear(c,numint);
						 break;
                case 84: control>>F84;
						 clear(c,numint);
						 break;
                case 85: control>>F85;
						 clear(c,numint);
						 break;
                case 150: control>>F150;
						 clear(c,numint);
						 break;
				case 151: control>>F151;
						 clear(c,numint);
						 break;
                 case 300: control>>F300;
						 clear(c,numint);
						 break;
				case 305: control>>F305;
						 clear(c,numint);
						 break;
				case 310: control>>F310;
						 clear(c,numint);
						 break;
				case 321: control>>F321;
						 clear(c,numint);
						 break;
				case 322: control>>F322;
						 clear(c,numint);
						 break;
				case 323: control>>F323;
						 clear(c,numint);
						 break;
				case 350: control>>F350;
						 clear(c,numint);
						 break;
				case 360: control>>F360;
						 clear(c,numint);
						 break;
				case 361: control>>F361;
						 clear(c,numint);
						 break;
				case 362: control>>F362;
						 clear(c,numint);
						 break;
			    case 369: ++F369;
						 clear(c,numint);
						 break;
				case 370: ++F370;
						 clear(c,numint);
						 break;
				case 371: ++F371;
						 clear(c,numint);
						 break;
                 case 374: ++F374;
						 clear(c,numint);
						 break;
                 case 375: ++F375;
						 clear(c,numint);
						 break;
                 case 378: ++F378;
						 clear(c,numint);
						 break;
                 case 379: ++F379;
						 clear(c,numint);
						 break;
				case 380: control>>F380;
						 clear(c,numint);
						 break;
				case 381: control>>F381;
						 clear(c,numint);
						 break;
				case 382: control>>F382;
						 clear(c,numint);
						 break;
				case 390: ++F390;
						 clear(c,numint);
						 break;
				case 391: ++F391;
						 clear(c,numint);
						 break;
                 case 394: ++F394;
						 clear(c,numint);
						 break;
                 case 395: ++F395;
						 clear(c,numint);
						 break;
                 case 398: ++F398;
						 clear(c,numint);
						 break;
                 case 399: ++F399;
						 clear(c,numint);
						 break;
				}
				break;

		case 'G':control>>numint;
				switch(numint)
				{
                 case 3: control>>G3;
						 clear(c,numint);
						 break;
				case 10: control>>G10;
						 clear(c,numint);
						 break;
				case 11: control>>G11;
						 clear(c,numint);
						 break;
				case 12: control>>G12;
						 clear(c,numint);
						 break;
				case 20: control>>G20;
						 clear(c,numint);
						 break;
				case 21: control>>G21;
						 clear(c,numint);
						 break;
				case 22: control>>G22;
						 clear(c,numint);
						 break;
				case 30: control>>G30;
						 clear(c,numint);
						 break;
                case 40: control>>G40;
						 clear(c,numint);
						 break;
				}
				break;

        case 'H': control>>numint;
				switch(numint)
				{
				 case 1: control>>H1;
						 clear(c,numint);
						 break;
                case 2: control>>H2;
						 clear(c,numint);
						 break;
                case 3: control>>H3;
						 clear(c,numint);
						 break;
                case 4: control>>H4_beta1>>H4_beta2;
                        H4=1;
						 clear(c,numint);
						 break;
                case 9: control>>H9;
						 clear(c,numint);
						 break;
                case 10: control>>H10;
						 clear(c,numint);
						 break;
                case 15: control>>H15;
						 clear(c,numint);
						 break;
                case 50: control>>H50_1>>H50_2;
						 clear(c,numint);
						 break;
                case 51: control>>H51;
						 clear(c,numint);
						 break;
                case 52: control>>H52;
						 clear(c,numint);
						 break;
                case 53: control>>H53;
						 clear(c,numint);
						 break;
                case 54: control>>H54;
						 clear(c,numint);
						 break;
                case 55: control>>H55;
						 clear(c,numint);
						 break;
                case 56: control>>H56;
						 clear(c,numint);
						 break;
                case 57: control>>H57_1>>H57_2>>H57_3>>H57_4;
						 clear(c,numint);
						 break;
                case 58: control>>H58_1>>H58_2>>H58_3>>H58_4;
						 clear(c,numint);
						 break;
                case 61: control>>H61_T;
                        H61=1;
						 clear(c,numint);
						 break;
                case 62: control>>H62_T;
                        H62=1;
						 clear(c,numint);
						 break;
                case 63: control>>H63_T;
                        H63=1;
						 clear(c,numint);
						 break;
                case 64: control>>H64_T;
                        H64=1;
						 clear(c,numint);
						 break;
                case 65: control>>H65_T;
                        H65=1;
						 clear(c,numint);
						 break;
                case 66: control>>H66_T;
                        H66=1;
						 clear(c,numint);
						 break;


				}
				break;

		case 'I': control>>numint;
				switch(numint)
				{
				case 10: control>>I10;
						 clear(c,numint);
						 break;
                case 11: control>>I11;
						 clear(c,numint);
						 break;
                case 12: control>>I12;
						 clear(c,numint);
						 break;
                case 13: control>>I13;
						 clear(c,numint);
						 break;
                case 21: control>>I21;
						 clear(c,numint);
						 break;
				 case 30: control>>I30;
						 clear(c,numint);
						 break;
				 case 40: control>>I40;
						 clear(c,numint);
						 break;
				 case 41: control>>I41;
						 clear(c,numint);
						 break;
                case 44: control>>I44;
						 clear(c,numint);
						 break;
                case 50: control>>I50;
						 clear(c,numint);
						 break;
                case 55: control>>I55;
						 clear(c,numint);
						 break;
                case 56: control>>I56;
						 clear(c,numint);
						 break;
                case 58: control>>I58_1>>I58_2;
						 clear(c,numint);
						 break;
                case 230: control>>I230;
						 clear(c,numint);
						 break;
                case 231: control>>I231;
						 clear(c,numint);
						 break;
                case 232: control>>I232;
						 clear(c,numint);
						 break;
                case 233: control>>I233;
						 clear(c,numint);
						 break;
                case 240: control>>I240;
						 clear(c,numint);
						 break;
                case 241: control>>I241;
						 clear(c,numint);
						 break;
				}
				break;

        case 'M': control>>numint;
				switch(numint)
				{
				case 10: control>>M10;
						 clear(c,numint);
						 break;

				}
				break;


		case 'N': control>>numint;
				switch(numint)
				{
				case 10: control>>N10;
						 clear(c,numint);
						 break;
                case 11: control>>N11;
						 clear(c,numint);
						 break;
				case 40: control>>N40;
						 clear(c,numint);
						 break;
				case 41: control>>N41;
						 clear(c,numint);
						 break;
				case 43: control>>N43;
						 clear(c,numint);
						 break;
				case 44: control>>N44;
						 clear(c,numint);
					     break;
				case 45: control>>N45;
						 clear(c,numint);
					     break;
				case 46: control>>N46;
						 clear(c,numint);
						 break;
				case 47: control>>N47;
						 clear(c,numint);
						 break;
                 case 48: control>>N48;
						 clear(c,numint);
						 break;
				case 49: control>>N49;
						 clear(c,numint);
						 break;
                 case 50: control>>N50;
						 clear(c,numint);
						 break;
                 case 60: control>>N60;
						 clear(c,numint);
						 break;
                 case 61: control>>N61;
						 clear(c,numint);
						 break;
				}
				break;

		case 'P': control>>numint;
				switch(numint)
				{
				case 10: control>>P10;
						 clear(c,numint);
						 break;
				case 11: control>>P11;
						 clear(c,numint);
						 break;
                case 12: control>>P12;
						 clear(c,numint);
						 break;
				case 14: control>>P14;
						 clear(c,numint);
						 break;
                case 15: control>>P15;
						 clear(c,numint);
						 break;
				case 18: control>>P18;
						 clear(c,numint);
						 break;
				case 20: control>>P20;
						 clear(c,numint);
						 break;
                case 21: control>>P21;
						 clear(c,numint);
						 break;
                case 22: control>>P22;
						 clear(c,numint);
						 break;
				case 23: control>>P23;
						 clear(c,numint);
						 break;
                case 24: control>>P24;
						 clear(c,numint);
						 break;
                case 25: control>>P25;
						 clear(c,numint);
						 break;
				case 26: control>>P26;
						 clear(c,numint);
						 break;
				case 27: control>>P27;
						 clear(c,numint);
						 break;
				case 28: control>>P28;
						 clear(c,numint);
						 break;
				case 29: control>>P29;
						 clear(c,numint);
						 break;
				case 30: control>>P30;
						 clear(c,numint);
						 break;
				case 34: control>>P34;
						 clear(c,numint);
						 break;
                case 35: ++P35;
						 clear(c,numint);
						 break;
                case 40: control>>P40;
						 clear(c,numint);
						 break;
                case 41: control>>P41;
						 clear(c,numint);
						 break;
				 case 42: control>>P42;
						 clear(c,numint);
						 break;
                case 43: control>>P43_xs>>P43_xe>>P43_ys>>P43_ye;
                        P43=1;
						 clear(c,numint);
						 break;
                case 44: control>>P44;
						 clear(c,numint);
						 break;
                case 45: control>>P45;
						 clear(c,numint);
						 break;
                case 50: ++P50;
						 clear(c,numint);
						 break;
				case 51: ++P51;
						 clear(c,numint);
						 break;
                case 52: ++P52;
						 clear(c,numint);
						 break;
				case 53: control>>P53;
						 clear(c,numint);
						 break;
				case 54: control>>P54;
						 clear(c,numint);
						 break;
				case 55: control>>P55;
						 clear(c,numint);
						 break;
				case 56: ++P56;
						 clear(c,numint);
						 break;
                case 57: ++P57;
						 clear(c,numint);
						 break;
                case 58: ++P58;
						 clear(c,numint);
						 break;
				case 59: control>>P59;
						 clear(c,numint);
						 break;
				case 61: ++P61;
						 clear(c,numint);
						 break;
				case 62: ++P62;
						 clear(c,numint);
						 break;
               case 63: ++P63;
						 clear(c,numint);
						 break;
                case 64: ++P64;
						 clear(c,numint);
						 break;
				case 66: control>>P66;
						 clear(c,numint);
						 break;
				case 67: ++P67;
						 clear(c,numint);
						 break;
                 case 68: ++P68;
						 clear(c,numint);
						 break;
                case 71: control>>P71;
						 clear(c,numint);
						 break;
                case 72: control>>P72;
						 clear(c,numint);
						 break;
                case 73: control>>P73;
						 clear(c,numint);
						 break;
                case 74: control>>P74;
						 clear(c,numint);
						 break;
                case 75: control>>P75;
						 clear(c,numint);
						 break;
                case 76: control>>P76;
						 clear(c,numint);
						 break;
                case 77: control>>P77;
						 clear(c,numint);
						 break;
                case 78: control>>P78;
						 clear(c,numint);
						 break;
			    case 79: control>>P79;
						 clear(c,numint);
						 break;
                case 81: ++P81;
						 clear(c,numint);
						 break;
               case 82: control>>P82;
						 clear(c,numint);
						 break;
				case 85: ++P85;
						 clear(c,numint);
						 break;
				case 91: control>>P91;
						 clear(c,numint);
						 break;
               case 92: control>>P92;
						 clear(c,numint);
						 break;
				case 101: control>>P101_xm>>P101_ym>>P101_zs>>P101_ze>>P101_r1>>P101_r2;
					      P101=1;
                          clear(c,numint);
						  break;
               case 120: control>>P120;
						 clear(c,numint);
						 break;
                case 110: control>>P110;
						 clear(c,numint);
						 break;
                case 111: control>>P111;
						 clear(c,numint);
						 break;
				case 121: ++P121;
						 clear(c,numint);
						 break;
				case 122: control>>P122;
						 clear(c,numint);
						 break;
			    case 123: ++P123;
						 clear(c,numint);
						 break;
				case 124: ++P124;
						 clear(c,numint);
						 break;
				case 125: ++P125;
						 clear(c,numint);
						 break;
				case 126: control>>P126;
						 clear(c,numint);
						 break;
				case 151: control>>P151;
						 clear(c,numint);
						 break;
				case 152: control>>P152;
						 clear(c,numint);
						 break;
				case 180: control>>P180;
						 clear(c,numint);
						 break;
				case 181: control>>P181;
						 clear(c,numint);
						 break;
				case 182: control>>P182;
						 clear(c,numint);
						 break;
               case 184: ++P184;
						 clear(c,numint);
						 break;
               case 185: ++P185;
						 clear(c,numint);
						 break;
                case 190: control>>P190;
						 clear(c,numint);
						 break;
				case 191: control>>P191;
						 clear(c,numint);
						 break;
				case 192: control>>P192;
						 clear(c,numint);
						 break;
               case 194: ++P194;
						 clear(c,numint);
						 break;
               case 195: ++P195;
						 clear(c,numint);
						 break;
               case 230: ++P230;
						 clear(c,numint);
						 break;
                case 240: ++P240;
						 clear(c,numint);
						 break;
				case 351: ++P351;
						 clear(c,numint);
						 break;
				case 352: ++P352;
						 clear(c,numint);
						 break;
				}
				break;
                
                
         case 'Q': control>>numint;
				switch(numint)
				{
				case 10: control>>Q10;
						 clear(c,numint);
						 break;
                 case 21: control>>Q21;
						 clear(c,numint);
						 break;
                 case 22: control>>Q22;
						 clear(c,numint);
						 break;
                 case 23: control>>Q23;
						 clear(c,numint);
						 break;
                 case 24: control>>Q24;
						 clear(c,numint);
						 break;
                 case 25: control>>Q25;
						 clear(c,numint);
						 break;
                 case 29: control>>Q29;
						 clear(c,numint);
						 break;
                 case 31: control>>Q31;
						 clear(c,numint);
						 break;
                 case 41: control>>Q41;
						 clear(c,numint);
						 break;
				 case 43: control>>Q43;
						 clear(c,numint);
						 break;
                 case 101: control>>Q101;
						 clear(c,numint);
						 break;
                 case 110: ++Q110;
						 clear(c,numint);
						 break;
                 case 111: control>>Q111_x;
                           Q111=1;
						 clear(c,numint);
						 break;
                 case 112: control>>Q112_y;
                           Q112=1;
						 clear(c,numint);
						 break;
                 case 113: control>>Q113_z;
                           Q113=1;
						 clear(c,numint);
						 break;
                 case 180: control>>Q180;
						 clear(c,numint);
						 break;
				case 181: control>>Q181;
						 clear(c,numint);
						 break;
				case 182: control>>Q182;
						 clear(c,numint);
						 break;

				}
				break;
                

        case 'S': control>>numint;
				switch(numint)
				{
				case 10: control>>S10;
						 clear(c,numint);
						 break;
                case 11: control>>S11;
						 clear(c,numint);
						 break;
                case 12: control>>S12;
						 clear(c,numint);
						 break;
                case 13: control>>S13;
						 clear(c,numint);
						 break;
                case 14: control>>S14;
						 clear(c,numint);
						 break;
                case 15: control>>S15;
						 clear(c,numint);
						 break;
                case 16: control>>S16;
						 clear(c,numint);
						 break;
               case 17: control>>S17;
						 clear(c,numint);
						 break;
				case 19: control>>S19;
						 clear(c,numint);
						 break;
				case 20: control>>S20;
						 clear(c,numint);
						 break;
				case 21: control>>S21;
						 clear(c,numint);
						 break;
                case 22: control>>S22;
						 clear(c,numint);
						 break;
                case 23: control>>S23_val;
                           S23=1;
						 clear(c,numint);
						 break;
                case 24: control>>S24;
						 clear(c,numint);
						 break;
                case 26: control>>S26_a>>S26_b;
						 clear(c,numint);
						 break;
                case 27: control>>S27;
						 clear(c,numint);
						 break;
                case 30: control>>S30;
						 clear(c,numint);
						 break;
                case 32: control>>S32;
						 clear(c,numint);
						 break;
                case 33: control>>S33;
						 clear(c,numint);
						 break;
                case 34: control>>S34;
						 clear(c,numint);
						 break;
				case 37: control>>S37;
						 clear(c,numint);
						 break;
                case 41: control>>S41;
						 clear(c,numint);
						 break;
			    case 42: control>>S42;
						 clear(c,numint);
						 break;
                case 43: control>>S43;
						 clear(c,numint);
						 break;
                case 44: control>>S44;
						 clear(c,numint);
						 break;
				 case 45: control>>S45;
						 clear(c,numint);
						 break;
				 case 46: control>>S46;
						 clear(c,numint);
						 break;
				 case 47: control>>S47;
						 clear(c,numint);
						 break;
				 case 48: control>>S48;
						 clear(c,numint);
						 break;
                case 50: control>>S50;
						 clear(c,numint);
						 break;
                case 57: control>>S57;
						 clear(c,numint);
						 break;
                case 60: control>>S60;
						 clear(c,numint);
						 break;
                case 71: control>>S71;
						 clear(c,numint);
						 break;
                case 72: control>>S72;
						 clear(c,numint);
						 break;
				 case 73: ++S73;
						 clear(c,numint);
						 break;
                case 77: control>>S77_xs>>S77_xe;
                          S77=1;
						 clear(c,numint);
						 break;
                case 78: control>>S78;
						 clear(c,numint);
						 break;
                case 79: control>>S79;
						 clear(c,numint);
						 break;
                case 80: control>>S80;
						 clear(c,numint);
						 break;
                case 81: control>>S81;
						 clear(c,numint);
						 break;
                case 82: control>>S82;
						 clear(c,numint);
						 break;
                case 83: control>>S83;
						 clear(c,numint);
						 break;
                case 84: control>>S84;
						 clear(c,numint);
						 break;
                case 90: control>>S90;
						 clear(c,numint);
						 break;
                case 91: control>>S91;
						 clear(c,numint);
						 break;
                case 92: control>>S92;
						 clear(c,numint);
						 break;
                case 93: control>>S93;
						 clear(c,numint);
						 break;
				case 100: control>>S100;
						 clear(c,numint);
						 break;
				case 101: control>>S101;
						 clear(c,numint);
						 break;
               case 116: control>>S116;
						 clear(c,numint);
						 break;
				}
				break;

		case 'T': control>>numint;
				switch(numint)
				{
				case 10: control>>T10;
						 clear(c,numint);
						 break;
                case 11: control>>T11;
						 clear(c,numint);
						 break;
                case 12: control>>T12;
						 clear(c,numint);
						 break;
                case 21: control>>T21;
						 clear(c,numint);
						 break;
                case 31: control>>T31;
						 clear(c,numint);
						 break;
				case 32: control>>T32;
						 clear(c,numint);
						 break;
                 case 33: control>>T33;
						 clear(c,numint);
						 break;
				case 35: control>>T35;
						 clear(c,numint);
						 break;
				case 36: control>>T36;
						 clear(c,numint);
						 break;
				case 37: control>>T37;
						 clear(c,numint);
						 break;
                case 38: control>>T38;
						 clear(c,numint);
						 break;
                case 39: control>>T39;
						 clear(c,numint);
						 break;
                case 41: control>>T41;
						 clear(c,numint);
						 break;
                case 42: control>>T42;
						 clear(c,numint);
						 break;
                case 43: control>>T43;
						 clear(c,numint);
						 break;
				}
				break;

		case 'W': control>>numint;
				switch(numint)
				{
				case  1: control>>W1;
						 clear(c,numint);
						 break;
				case  2: control>>W2;
						 clear(c,numint);
						 break;
				case  3: control>>W3;
						 clear(c,numint);
						 break;
				case  4: control>>W4;
						 clear(c,numint);
						 break;
				case  5: control>>W5;
						 clear(c,numint);
						 break;
                case  6: control>>W6;
						 clear(c,numint);
						 break;
                case  7: control>>W7;
						 clear(c,numint);
						 break;
				case 10: control>>W10;
						 clear(c,numint);
						 break;
               case 11: control>>W11_u>>W11_v>>W11_w;
                        W11=1;
						 clear(c,numint);
						 break;
               case 12: control>>W12_u>>W12_v>>W12_w;
                        W12=1;
						 clear(c,numint);
						 break;
               case 13: control>>W13_u>>W13_v>>W13_w;
                        W13=1;
						 clear(c,numint);
						 break;
               case 14: control>>W14_u>>W14_v>>W14_w;
                        W14=1;
						 clear(c,numint);
						 break;
               case 15: control>>W15_u>>W15_v>>W15_w;
                        W15=1;
						 clear(c,numint);
						 break;
               case 16: control>>W16_u>>W16_v>>W16_w;
                        W16=1;
						 clear(c,numint);
						 break;
				case 20: control>>W20;
						 clear(c,numint);
						 break;
				case 21: control>>W21;
						 clear(c,numint);
						 break;
				case 22: control>>W22;
						 clear(c,numint);
						 break;
               case 29: control>>W29_x>>W29_y>>W29_z;
						 clear(c,numint);
						 break;
				case 30: control>>W30;
						 clear(c,numint);
						 break;
				case 31: control>>W31;
						 clear(c,numint);
						 break;
               case 41: ++W41;
						 clear(c,numint);
						 break;
                case 50: control>>W50;
                        W50_air=1;
						 clear(c,numint);
						 break;
                case 90: control>>W90;
						 clear(c,numint);
						 break;
                case 95: control>>W95;
						 clear(c,numint);
						 break;
                case 96: control>>W96;
						 clear(c,numint);
						 break;
                case 97: control>>W97;
						 clear(c,numint);
						 break;
                case 98: control>>W98;
						 clear(c,numint);
						 break;
                case 101: control>>W101;
						 clear(c,numint);
						 break;
                case 102: control>>W102_phi>>W102_c;
						 clear(c,numint);
						 break;
                case 103: control>>W103;
						 clear(c,numint);
						 break;
                case 104: control>>W104;
						 clear(c,numint);
						 break;
                case 110: control>>W110;
						 clear(c,numint);
						 break;
                case 111: control>>W111;
						 clear(c,numint);
						 break;
                case 112: control>>W112;
						 clear(c,numint);
						 break;
				}
				break;

		case 'X': control>>numint;
				switch(numint)
				{
				case  10: control>>X10;
						 clear(c,numint);
						 break;
				case  11: control>>X11_u>>X11_v>>X11_w>>X11_p>>X11_q>>X11_r;
						 clear(c,numint);
						 break;
				case  12: control>>X12;
						 clear(c,numint);
						 break;
                 case  14: control>>X14;
						 clear(c,numint);
						 break;
                 case  15: control>>X15;
						 clear(c,numint);
						 break;
				case  19: control>>X19;
						 clear(c,numint);
						 break;
				case  21: control>>X21_d;
						 X21=1;
						 clear(c,numint);
						 break;
				case  22: control>>X22_m;
						 X22=1;
						 clear(c,numint);
						 break;
				case  23: control>>X23_x>>X23_y>>X23_z;
						 X23=1;
						 clear(c,numint);
						 break;
				case  24: control>>X24_Ix>>X24_Iy>>X24_Iz;
						 X24=1;
						 clear(c,numint);
						 break;
				case  25: control>>X25_Cp>>X25_Cq>>X25_Cr;
						 clear(c,numint);
						 break;
                case  26: control>>X26_Cu>>X26_Cv>>X26_Cw;
						 clear(c,numint);
						 break;
			    case  31: control>>X31;
						 clear(c,numint);
						 break;
				case  32: control>>X32;
						 clear(c,numint);
						 break;
				case  33: control>>X33;
						 clear(c,numint);
						 break;
                case  34: control>>X34;
						 clear(c,numint);
						 break;
                case  39: control>>X39;
						 clear(c,numint);
						 break;
				case  40: control>>X40;
						 clear(c,numint);
						 break;
               case  41: control>>X41;
						 clear(c,numint);
						 break;
               case  42: control>>X42;
						 clear(c,numint);
						 break;
               case  43: control>>X43;
						 clear(c,numint);
						 break;
               case  44: control>>X44;
						 clear(c,numint);
						 break;
                case  45: control>>X45;
						 clear(c,numint);
						 break;
                case  46: control>>X46;
						 clear(c,numint);
						 break;
                case  47: control>>X47;
						 clear(c,numint);
						 break;
                case  48: control>>X48;
						 clear(c,numint);
						 break;
                case  49: control>>X49;
						 clear(c,numint);
						 break;
               case  50: control>>X50;
						 clear(c,numint);
						 break;
				case  100: control>>X100_x>>X100_y>>X100_z;
						 X100=1;
						 clear(c,numint);
						 break;
				case  101: control>>X101_phi>>X101_theta>>X101_psi;
						 X101=1;
						 clear(c,numint);
						 break;
				case  102: control>>X102_u>>X102_v>>X102_w;
						 X102=1;
						 clear(c,numint);
						 break;
				case  103: control>>X103_p>>X103_q>>X103_r;
						 X103=1;
						 clear(c,numint);
						 break;
				case  110: ++X110;
						 clear(c,numint);
						 break;
				case  120: control>>X120_rad>>X120_xc>>X120_yc>>X120_zc;
						 X120=1;
						 clear(c,numint);
						 break;
				case  131: control>>X131_rad>>X131_h>>X131_xc>>X131_yc>>X131_zc;
						 X131=1;
						 clear(c,numint);
						 break;
				case  132: control>>X132_rad>>X132_h>>X132_xc>>X132_yc>>X132_zc;
						 X132=1;
						 clear(c,numint);
						 break;
				case  133: control>>X133_rad>>X133_h>>X133_xc>>X133_yc>>X133_zc;
						 X133=1;
						 clear(c,numint);
						 break;
				case  153: control>>X153_xs>>X153_xe>>X153_ys>>X153_ye>>X153_zs>>X153_ze;
						 X153=1;
						 clear(c,numint);
						 break;
                case  163: ++X163;
						 clear(c,numint);
						 break;
                case  164: ++X164;
						 clear(c,numint);
						 break;
                case  180: control>>X180;
						 clear(c,numint);
						 break;
                case  181: control>>X181_x>>X181_y>>X181_z;
                         X181=1;
						 clear(c,numint);
						 break;
                case  182: control>>X182_x>>X182_y>>X182_z;
                         X182=1;
						 clear(c,numint);
						 break;
                case  183: control>>X183_x>>X183_y>>X183_z>>X183_phi>>X183_theta>>X183_psi;
                         X183=1;
						 clear(c,numint);
						 break;
                case  184: control>>X184;
						 clear(c,numint);
						 break;
                case  205: control>>X205;
						 clear(c,numint);
						 break;
                case  206: control>>X206_ts>>X206_te;
                            X206=1;
						 clear(c,numint);
						 break;
                case  207: control>>X207_ts>>X207_te;
                            X207=1;
						 clear(c,numint);
						 break;
				case  210: control>>X210_u>>X210_v>>X210_w;
						 X210=1;
						 clear(c,numint);
						 break;
				case  211: control>>X211_p>>X211_q>>X211_r;
						 X211=1;
						 clear(c,numint);
						 break;
                case  221: control>>X221_xs>>X221_xe>>X221_ys>>X221_ye>>X221_zs>>X221_ze;
						 X221=1;
						 clear(c,numint);
						 break;
				case  310: control>>X310;
						 clear(c,numint);
						 break;
                case  311: ++X311;
						 clear(c,numint);
						 break;
                case  312: ++X312;
						 clear(c,numint);
						 break;
				case  313: control>>X313;
						 clear(c,numint);
						 break;
                case  314: X314 = 1; 
						 clear(c,numint);
						 break;
                case  315:  X315 = 1;
						 clear(c,numint);
						 break;
				case  320: ++X320;
                         B269=3;
						 clear(c,numint);
						 break;
				case  321: ++X321;
						 clear(c,numint);
						 break;
				case  323: control>>X323_m>>X323_d>>X323_l;
                         clear(c,numint);
						 break;
				case  324: ++X324;
						 clear(c,numint);
						 break;
				case  325: control>>X325_dt>>X325_relX>>X325_relY>>X325_relZ;
                         clear(c,numint);
						 break;
				case  400: control>>X400;
                         clear(c,numint);
						 break;
				case  401: control>>X401_p0>>X401_cl>>X401_cb>>X401_a;
                         clear(c,numint);
						 break;
				}
				break;

		case 'Y': control>>numint;
				switch(numint)
				{
               case 1: control>>Y1;
                        clear(c,numint);
                        break;
               case 2: control>>Y2;
                        clear(c,numint);
                        break;
               case 3: control>>Y3;
                        clear(c,numint);
                        break;
               case 4: control>>Y4;
                        clear(c,numint);
                        break;
               case 5: control>>Y5;
                        clear(c,numint);
                        break;
               case 40: control>>Y40;
						 clear(c,numint);
						 break;
				case 50: control>>Y50;
						 clear(c,numint);
						 break;
				case 60: control>>Y60;
						 clear(c,numint);
						 break;
               case 71: control>>Y71;
						 clear(c,numint);
						 break;
               case 72: control>>Y72;
						 clear(c,numint);
						 break;
               case 73: control>>Y73;
						 clear(c,numint);
						 break;
               case 74: control>>Y74;
						 clear(c,numint);
						 break;

				}
				break;

		case 'Z': control>>numint;
				switch(numint)
				{
				case  10: control>>Z10;
						 clear(c,numint);
						 break;
                case  11: ++Z11;
						 clear(c,numint);
						 break;
                }
		}

        ++count;

	}
        if(count>1e7)
        {
        cout<<endl;
        cout<<"!!! missing input parameter in ctrl.txt !!!"<<endl<<endl;
        cout<<"!!! please check the REEF3D User Guide !!!"<<endl<<endl<<endl<<endl;

        exit(0);
        }
	}
	control.close();
	control.clear();


    // re-read

	// B
	Darray(B70_val,B70);
	Darray(B70_dist,B70);
	Darray(B70_b,B70);
	Darray(B70_x,B70);
	Darray(B70_y,B70);

	Darray(B71_val,B71);
	Darray(B71_dist,B71);
	Darray(B71_b,B71);
	Darray(B71_x,B71);
	Darray(B71_y,B71);

	Darray(B106_b,B106);
	Darray(B106_x,B106);
	Darray(B106_y,B106);

	Darray(B107_xs,B107);
    Darray(B107_xe,B107);
	Darray(B107_ys,B107);
    Darray(B107_ye,B107);
    Darray(B107_d,B107);

    Darray(B108_xs,B108);
    Darray(B108_xe,B108);
	Darray(B108_ys,B108);
    Darray(B108_ye,B108);
    Darray(B108_d,B108);


	Darray(B240_C,B240);
	Darray(B240_D,B240);
	Darray(B240_xs,B240);
	Darray(B240_xe,B240);
	Darray(B240_ys,B240);
	Darray(B240_ye,B240);
	Darray(B240_zs,B240);
	Darray(B240_ze,B240);

    Darray(B270_xs,B270);
	Darray(B270_xe,B270);
	Darray(B270_ys,B270);
	Darray(B270_ye,B270);
	Darray(B270_zs,B270);
	Darray(B270_ze,B270);
    Darray(B270_n,B270);
    Darray(B270_d50,B270);
	Darray(B270_alpha,B270);
	Darray(B270_beta,B270);

	Darray(B274_xc,B274);
	Darray(B274_yc,B274);
	Darray(B274_zs,B274);
	Darray(B274_ze,B274);
	Darray(B274_r,B274);
    Darray(B274_n,B274);
    Darray(B274_d50,B274);
	Darray(B274_alpha,B274);
	Darray(B274_beta,B274);

    Darray(B281_xs,B281);
	Darray(B281_xe,B281);
	Darray(B281_ys,B281);
	Darray(B281_ye,B281);
	Darray(B281_zs,B281);
	Darray(B281_ze,B281);
    Darray(B281_n,B281);
    Darray(B281_d50,B281);
	Darray(B281_alpha,B281);
	Darray(B281_beta,B281);
    
    Darray(B282_xs,B282);
	Darray(B282_xe,B282);
	Darray(B282_ys,B282);
	Darray(B282_ye,B282);
	Darray(B282_zs,B282);
	Darray(B282_ze,B282);
    Darray(B282_n,B282);
    Darray(B282_d50,B282);
	Darray(B282_alpha,B282);
	Darray(B282_beta,B282);

    Darray(B291_xs,B291);
	Darray(B291_xe,B291);
	Darray(B291_ys,B291);
	Darray(B291_ye,B291);
	Darray(B291_zs,B291);
	Darray(B291_ze,B291);
    Darray(B291_d,B291);
    Darray(B291_n,B291);
    Darray(B291_d50,B291);
	Darray(B291_alpha,B291);
	Darray(B291_beta,B291);

    Darray(B310_xs,B310);
	Darray(B310_xe,B310);
	Darray(B310_ys,B310);
	Darray(B310_ye,B310);
	Darray(B310_zs,B310);
	Darray(B310_ze,B310);
    Darray(B310_N,B310);
    Darray(B310_D,B310);
	Darray(B310_Cd,B310);
    
    Darray(B321_xs,B321);
	Darray(B321_xe,B321);
	Darray(B321_ys,B321);
	Darray(B321_ye,B321);
	Darray(B321_zs,B321);
	Darray(B321_ze,B321);
    Darray(B321_N,B321);
    Darray(B321_D,B321);
	Darray(B321_Cd,B321);
    
    Darray(B322_xs,B322);
	Darray(B322_xe,B322);
	Darray(B322_ys,B322);
	Darray(B322_ye,B322);
	Darray(B322_zs,B322);
	Darray(B322_ze,B322);
    Darray(B322_N,B322);
    Darray(B322_D,B322);
	Darray(B322_Cd,B322);


    Iarray(B411_ID,B411);
    Darray(B411_Q,B411);

    Iarray(B412_ID,B412);
    Darray(B412_pressBC,B412);

    Iarray(B413_ID,B413);
    Darray(B413_h,B413);

    Iarray(B414_ID,B414);
    Darray(B414_Uio,B414);

    Iarray(B415_ID,B415);
    Darray(B415_U,B415);
    Darray(B415_V,B415);
    Darray(B415_W,B415);
    Iarray(B416_ID,B416);
    Darray(B416_alpha,B416);
    Iarray(B417_ID,B417);
    Darray(B417_Nx,B417);
    Darray(B417_Ny,B417);
    Darray(B417_Nz,B417);
    Iarray(B418_ID,B418);
    Iarray(B418_pio,B418);
    Iarray(B421_ID,B421);
    Iarray(B421_Q,B421);
    Iarray(B422_ID,B422);
    Iarray(B422_FSF,B422);


    Iarray(B440_ID,B440);
    Iarray(B440_face,B440);
    Darray(B440_xs,B440);
    Darray(B440_xe,B440);
    Darray(B440_ys,B440);
    Darray(B440_ye,B440);

    Iarray(B441_ID,B441);
    Iarray(B441_face,B441);
    Darray(B441_xs,B441);
    Darray(B441_xe,B441);
    Darray(B441_ys,B441);
    Darray(B441_ye,B441);
    Darray(B441_zs,B441);
    Darray(B441_ze,B441);

    Iarray(B442_ID,B442);
    Iarray(B442_face,B442);
    Darray(B442_xm,B442);
    Darray(B442_ym,B442);
    Darray(B442_zm,B442);
    Darray(B442_r,B442);

	// C
	Darray(C75_x,C75);
	Darray(C75_z,C75);
	Darray(C75_a,C75);
	Darray(C75_s,C75);
	Darray(C75_l,C75);
	Darray(C75_v,C75);


    // F
	Darray(F70_xs,F70);
	Darray(F70_xe,F70);

	Darray(F70_ys,F70);
	Darray(F70_ye,F70);

	Darray(F70_zs,F70);
	Darray(F70_ze,F70);


	Darray(F71_xs,F71);
	Darray(F71_xe,F71);

	Darray(F71_ys,F71);
	Darray(F71_ye,F71);

	Darray(F71_zs,F71);
	Darray(F71_ze,F71);


	Darray(F72_xs,F72);
	Darray(F72_xe,F72);

	Darray(F72_ys,F72);
	Darray(F72_ye,F72);

	Darray(F72_h,F72);
    
    if(F369>0)
	{
	Darray(F369_x,F369);   
	Darray(F369_z,F369);  
	
	Darray(F369_a,F369);  	
	Darray(F369_s,F369);  
	Darray(F369_l,F369);  
	Darray(F369_v,F369);  
	}
    
	if(F370>0)
	{
	Darray(F370_xs,F370);  
	Darray(F370_xe,F370);  
	
	Darray(F370_ys,F370);  
	Darray(F370_ye,F370);  
	
	Darray(F370_zs,F370);  
	Darray(F370_ze,F370);  
	}
	
	if(F371>0)
	{
	Darray(F371_xs,F371);  
	Darray(F371_xe,F371);  
	
	Darray(F371_ys,F371);  
	Darray(F371_ye,F371);  
	
	Darray(F371_zs,F371);  
	Darray(F371_ze,F371);  
	}
	
    if(F374>0)
	{
	Darray(F374_xc,F374);  
	Darray(F374_zc,F374);
    Darray(F374_r,F374);  
	}
    
    if(F375>0)
	{
	Darray(F375_xc,F375);  
	Darray(F375_zc,F375);
    Darray(F375_r,F375);  
	}
    
    if(F378>0)
	{
	Darray(F378_xc,F378);  
    Darray(F378_yc,F378);  
	Darray(F378_zc,F378);
    Darray(F378_r,F378);  
	}
	
    if(F379>0)
	{
	Darray(F379_xc,F379);  
    Darray(F379_yc,F379);  
	Darray(F379_zc,F379);
    Darray(F379_r,F379);  
	}
	
	if(F390>0)
	{
	Darray(F390_xs,F390);  
	Darray(F390_xe,F390);  
	
	Darray(F390_ys,F390);  
	Darray(F390_ye,F390);  
	
	Darray(F390_zs,F390);  
	Darray(F390_ze,F390);  
	}
	
	if(F391>0)
	{
	Darray(F391_xs,F391);  
	Darray(F391_xe,F391);  
	
	Darray(F391_ys,F391);  
	Darray(F391_ye,F391);  
	
	Darray(F391_zs,F391);  
	Darray(F391_ze,F391);  
	}
    
    if(F394>0)
	{
	Darray(F394_xc,F394);  
	Darray(F394_zc,F394);
    Darray(F394_r,F394);  
	}
    
    if(F395>0)
	{
	Darray(F395_xc,F395);  
	Darray(F395_zc,F395);
    Darray(F395_r,F395);  
	}
    
    if(F398>0)
	{
	Darray(F398_xc,F398);  
    Darray(F398_yc,F398);  
	Darray(F398_zc,F398);
    Darray(F398_r,F398);  
	}
	
    if(F399>0)
	{
	Darray(F399_xc,F399);  
    Darray(F399_yc,F399);  
	Darray(F399_zc,F399);
    Darray(F399_r,F399);  
	}

    // P
	Darray(P35_ts,P35);
	Darray(P35_te,P35);
	Darray(P35_dt,P35);

	Darray(P50_x,P50);
	Darray(P50_y,P50);

	Darray(P51_x,P51);
	Darray(P51_y,P51);

	Darray(P52_y,P52);
	Darray(P56_x,P56);

	Darray(P58_x,P58);
	Darray(P58_y,P58);
	Darray(P58_T,P58);
    
    Darray(P61_x,P61);
	Darray(P61_y,P61);
	Darray(P61_z,P61);

	Darray(P62_xs,P62);
	Darray(P62_ys,P62);
	Darray(P62_zs,P62);
	Darray(P62_xe,P62);
	Darray(P62_ye,P62);
	Darray(P62_ze,P62);

    Darray(P63_x,P63);
	Darray(P63_y,P63);
    
    Darray(P64_x,P64);
	Darray(P64_y,P64);
	Darray(P64_z,P64);

	Darray(P67_x,P67);
    
    Darray(P68_x,P68);
    Darray(P68_zs,P68);
    Darray(P68_ze,P68);

	Darray(P81_xs,P81);
	Darray(P81_ys,P81);
	Darray(P81_zs,P81);
	Darray(P81_xe,P81);
	Darray(P81_ye,P81);
	Darray(P81_ze,P81);

	Darray(P85_x,P85);
	Darray(P85_y,P85);
	Darray(P85_r,P85);
	Darray(P85_cd,P85);
	Darray(P85_cm,P85);

	Darray(P121_x,P121);
	Darray(P121_y,P121);

	Darray(P123_y,P123);

	Darray(P124_x,P124);

	Darray(P125_x,P125);
	Darray(P125_y,P125);

    Iarray(P184_its,P184);
	Iarray(P184_ite,P184);
	Iarray(P184_dit,P184);

    Darray(P185_ts,P185);
	Darray(P185_te,P185);
	Darray(P185_dt,P185);
    
    Iarray(P194_its,P194);
	Iarray(P194_ite,P194);
	Iarray(P194_dit,P194);

    Darray(P195_ts,P195);
	Darray(P195_te,P195);
	Darray(P195_dt,P195);

    Darray(P230_x,P230);
    Darray(P240_x,P240);

	Darray(P351_x,P351);
	Darray(P351_y,P351);
	Darray(P352_x,P352);
	Darray(P352_y,P352);
    
    // Q
	Darray(Q110_xs,Q110);
	Darray(Q110_ys,Q110);
	Darray(Q110_zs,Q110);
	Darray(Q110_xe,Q110);
	Darray(Q110_ye,Q110);
	Darray(Q110_ze,Q110);

	// S
	Darray(S73_val,S73);
	Darray(S73_dist,S73);
	Darray(S73_b,S73);
	Darray(S73_x,S73);
	Darray(S73_y,S73);

    // W
    Darray(W41_xc,W41);
    Darray(W41_yc,W41);
    Darray(W41_zs,W41);
    Darray(W41_ze,W41);
    Darray(W41_vel,W41);
    Darray(W41_beta,W41);

	// X
	Darray(X110_xs,X110);
	Darray(X110_ys,X110);
	Darray(X110_zs,X110);
	Darray(X110_xe,X110);
	Darray(X110_ye,X110);
	Darray(X110_ze,X110);

    Darray(X163_x1,X163);
    Darray(X163_y1,X163);
    Darray(X163_z1,X163);
    Darray(X163_x2,X163);
    Darray(X163_y2,X163);
    Darray(X163_z2,X163);
    Darray(X163_x3,X163);
    Darray(X163_y3,X163);
    Darray(X163_z3,X163);
    Darray(X163_x4,X163);
    Darray(X163_y4,X163);
    Darray(X163_z4,X163);
    Darray(X163_x5,X163);
    Darray(X163_y5,X163);
    Darray(X163_z5,X163);
    Darray(X163_x6,X163);
    Darray(X163_y6,X163);
    Darray(X163_z6,X163);

    Darray(X164_x1,X164);
    Darray(X164_y1,X164);
    Darray(X164_z1,X164);
    Darray(X164_x2,X164);
    Darray(X164_y2,X164);
    Darray(X164_z2,X164);
    Darray(X164_x3,X164);
    Darray(X164_y3,X164);
    Darray(X164_z3,X164);
    Darray(X164_x4,X164);
    Darray(X164_y4,X164);
    Darray(X164_z4,X164);
    Darray(X164_x5,X164);
    Darray(X164_y5,X164);
    Darray(X164_z5,X164);
    Darray(X164_x6,X164);
    Darray(X164_y6,X164);
    Darray(X164_z6,X164);
    Darray(X164_x7,X164);
    Darray(X164_y7,X164);
    Darray(X164_z7,X164);
    Darray(X164_x8,X164);
    Darray(X164_y8,X164);
    Darray(X164_z8,X164);

    if (X311 > 0)
    {
        Darray(X311_xs,X311);
        Darray(X311_ys,X311);
        Darray(X311_zs,X311);
        Darray(X311_xe,X311);
        Darray(X311_ye,X311);
        Darray(X311_ze,X311);
        Darray(X311_w,X311);
        Darray(X311_rho_c,X311);
        Darray(X311_EA,X311);
        Darray(X311_d,X311);
        Darray(X311_l,X311);
        Darray(X311_H,X311);
        Darray(X311_P,X311);
        Darray(X311_facT,X311);

        mooring_count = X311;

        Darray(X314_T,X311);
        Darray(X315_t,X311);
    }
    else
    {
        Darray(X311_xs,X312);
        Darray(X311_ys,X312);
        Darray(X311_zs,X312);
        Darray(X311_xe,X312);
        Darray(X311_ye,X312);
        Darray(X311_ze,X312);
        Darray(X312_k,X312);
        Darray(X312_T0,X312);

        Darray(X311_w,X311);
        Darray(X311_rho_c,X311);
        Darray(X311_EA,X311);
        Darray(X311_d,X311);
        Darray(X311_l,X311);
        Darray(X311_H,X311);
        Darray(X311_P,X311);
        Darray(X311_facT,X311);

        mooring_count = X312;

        Darray(X314_T,X312);
        Darray(X315_t,X312);
    }

    if (X321 > 0)
    {
		Iarray(X320_type,X320);

        Darray(X321_Sn,X321);
		Darray(X321_d,X321);
        Darray(X321_lambda,X321);
		Darray(X321_dk,X321);
		Darray(X321_rho,X321);
		Darray(X321_nd,X321);
		Darray(X321_nl,X321);

        Darray(X322_D,X321);
        Darray(X322_L,X321);
        Darray(X322_x0,X321);
        Darray(X322_y0,X321);
        Darray(X322_z0,X321);
        Darray(X322_phi,X321);
        Darray(X322_theta,X321);
        Darray(X322_psi,X321);

        net_count = X321;
    }

    Darray(X324_x,X324);
	Darray(X324_y,X324);
	Darray(X324_z,X324);

    if (Z11 > 0)
    {
        Darray(Z11_x,Z11);
        Darray(Z11_y,Z11);
        Darray(Z11_z,Z11);
        Darray(Z11_l,Z11);
        Darray(Z11_w,Z11);
        Darray(Z11_t,Z11);
        Darray(Z11_rho,Z11);
        Darray(Z11_e,Z11);
        Darray(Z11_ix,Z11);
        Darray(Z11_iy,Z11);
        Darray(Z11_iz,Z11);
        Darray(Z11_nu,Z11);
        Darray(Z11_n,Z11);

        FSI_count = Z11;
    }
    
	int countB70=0;
	int countB71=0;
	int countB106=0;
	int countB107=0;
    int countB108=0;
    int countB231=0;
    int countB232=0;
	int countB240=0;
    int countB270=0;
    int countB274=0;
    int countB281=0;
    int countB282=0;
    int countB291=0;
    int countB310=0;
    int countB321=0;
    int countB322=0;
    int countB411=0;
    int countB412=0;
    int countB413=0;
    int countB414=0;
    int countB415=0;
    int countB416=0;
    int countB417=0;
    int countB418=0;
    int countB421=0;
    int countB422=0;
    int countB440=0;
    int countB441=0;
    int countB442=0;
    int countC75=0;
	int countF70=0;
	int countF71=0;
	int countF72=0;
    int countF369=0;
	int countF370=0;
	int countF371=0;
	int countF374=0;
    int countF375=0;
    int countF378=0;
    int countF379=0;
	int countF390=0;
	int countF391=0;
    int countF394=0;
    int countF395=0;
    int countF398=0;
    int countF399=0;
	int countP35=0;
    int countP50=0;
	int countP51=0;
    int countP52=0;
	int countP56=0;
    int countP58=0;
	int countP61=0;
	int countP62=0;
    int countP63=0;
    int countP64=0;
	int countP67=0;
    int countP68=0;
	int countP81=0;
	int countP85=0;
	int countP121=0;
	int countP123=0;
	int countP124=0;
	int countP125=0;
    int countP184=0;
    int countP185=0;
    int countP194=0;
    int countP195=0;
    int countP230=0;
    int countP240=0;
	int countP351=0;
	int countP352=0;
    int countQ110=0;
	int countS73=0;
    int countW41=0;
	int countX110=0;
    int countX163=0;
    int countX164=0;
    int countX311=0;
    int countX312=0;
    int countX320=0;
    int countX321=0;
    int countX322=0;
	int countX324=0;
	int countZ11=0;

	control.open("ctrl.txt", ios_base::in);
	while(!control.eof())
	{
		control>>c;
		if (c == '/')
			control.ignore(1000, '\n');
		else
		{
			switch(c)
			{
			case 'B': control>>numint;
				switch(numint)
				{
				case 70: control>>B70_val[countB70]>>B70_dist[countB70]>>B70_b[countB70]>>B70_x[countB70]>>B70_y[countB70];
                        ++countB70;
						 clear(c,numint);
						 break;
				case 71: control>>B71_val[countB71]>>B71_dist[countB71]>>B71_b[countB71]>>B71_x[countB71]>>B71_y[countB71];
                        ++countB71;
						 clear(c,numint);
						 break;
				case 106: control>>B106_b[countB106]>>B106_x[countB106]>>B106_y[countB106];
                        ++countB106;
						 clear(c,numint);
						 break;
				case 107: control>>B107_xs[countB107]>>B107_xe[countB107]>>B107_ys[countB107]>>B107_ye[countB107]>>B107_d[countB107];
                        ++countB107;
						 clear(c,numint);
						 break;
                case 108: control>>B108_xs[countB108]>>B108_xe[countB108]>>B108_ys[countB108]>>B108_ye[countB108]>>B108_d[countB108];
                        ++countB108;
						 clear(c,numint);
						 break;
                case 240: control>>B240_C[countB240]>>B240_D[countB240]>>B240_xs[countB240]>>B240_xe[countB240]>>B240_ys[countB240]>>B240_ye[countB240]>>B240_zs[countB240]>>B240_ze[countB240];
                        ++countB240;
						 clear(c,numint);
						 break;
                case 270: control>>B270_xs[countB270]>>B270_xe[countB270]>>B270_ys[countB270]>>B270_ye[countB270]>>B270_zs[countB270]>>B270_ze[countB270]>>B270_n[countB270]>>B270_d50[countB270]>>B270_alpha[countB270]>>B270_beta[countB270];
                        ++countB270;
						 clear(c,numint);
						 break;
                case 274: control>>B274_xc[countB274]>>B274_yc[countB274]>>B274_zs[countB274]>>B274_ze[countB274]>>B274_r[countB274]>>B274_n[countB274]>>B274_d50[countB274]>>B274_alpha[countB274]>>B274_beta[countB274];
                        ++countB274;
						 clear(c,numint);
						 break;
                case 281: control>>B281_xs[countB281]>>B281_xe[countB281]>>B281_ys[countB281]>>B281_ye[countB281]>>B281_zs[countB281]>>B281_ze[countB281]>>B281_n[countB281]>>B281_d50[countB281]>>B281_alpha[countB281]>>B281_beta[countB281];
                        ++countB281;
						 clear(c,numint);
						 break;
                case 282: control>>B282_xs[countB282]>>B282_xe[countB282]>>B282_ys[countB282]>>B282_ye[countB282]>>B282_zs[countB282]>>B282_ze[countB282]>>B282_n[countB282]>>B282_d50[countB282]>>B282_alpha[countB282]>>B282_beta[countB282];
                        ++countB282;
						 clear(c,numint);
						 break;
                case 291: control>>B291_xs[countB291]>>B291_xe[countB291]>>B291_ys[countB291]>>B291_ye[countB291]>>B291_zs[countB291]>>B291_ze[countB291]>>B291_d[countB291]>>B291_n[countB291]>>B291_d50[countB291]>>B291_alpha[countB291]>>B291_beta[countB291];
                        ++countB291;
						 clear(c,numint);
						 break;
                case 310: control>>B310_xs[countB310]>>B310_xe[countB310]>>B310_ys[countB310]>>B310_ye[countB310]>>B310_zs[countB310]>>B310_ze[countB310]>>B310_N[countB310]>>B310_D[countB310]>>B310_Cd[countB310];
                        ++countB310;
						 clear(c,numint);
						 break;
                case 321: control>>B321_xs[countB321]>>B321_xe[countB321]>>B321_ys[countB321]>>B321_ye[countB321]>>B321_zs[countB321]>>B321_ze[countB321]>>B321_N[countB321]>>B321_D[countB321]>>B321_Cd[countB321];
                        ++countB321;
						 clear(c,numint);
						 break;
                case 322: control>>B322_xs[countB322]>>B322_xe[countB322]>>B322_ys[countB322]>>B322_ye[countB322]>>B322_zs[countB322]>>B322_ze[countB322]>>B322_N[countB322]>>B322_D[countB322]>>B322_Cd[countB322];
                        ++countB322;
						 clear(c,numint);
						 break;
                case 411: control>>B411_ID[countB411]>>B411_Q[countB411];
                        ++countB411;
						 clear(c,numint);
						 break;
                case 412: control>>B412_ID[countB412]>>B412_pressBC[countB412];
                        ++countB412;
						 clear(c,numint);
						 break;
                case 413: control>>B413_ID[countB413]>>B413_h[countB413];
                        ++countB413;
						 clear(c,numint);
						 break;
                case 414: control>>B414_ID[countB414]>>B414_Uio[countB414];
                        ++countB414;
						 clear(c,numint);
						 break;
                case 415: control>>B415_ID[countB415]>>B415_U[countB415]>>B415_V[countB415]>>B415_W[countB415];
                        ++countB415;
						 clear(c,numint);
						 break;
                case 416: control>>B416_ID[countB416]>>B416_alpha[countB416];
                        ++countB416;
						 clear(c,numint);
						 break;
                case 417: control>>B417_ID[countB417]>>B417_Nx[countB417]>>B417_Ny[countB417]>>B417_Nz[countB417];
                        ++countB417;
						 clear(c,numint);
						 break;
                case 418: control>>B418_ID[countB418]>>B418_pio[countB418];
                        ++countB418;
						 clear(c,numint);
						 break;
                case 421: control>>B421_ID[countB421]>>B421_Q[countB421];
                        ++countB421;
						 clear(c,numint);
						 break;
                case 422: control>>B422_ID[countB422]>>B422_FSF[countB422];
                        ++countB422;
						 clear(c,numint);
						 break;
                case 440: control>>B440_ID[countB440]>>B440_face[countB440]>>B440_xs[countB440]>>B440_xe[countB440]>>B440_ys[countB440]>>B440_ye[countB440];
                        ++countB440;
						 clear(c,numint);
						 break;
                case 441: control>>B441_ID[countB441]>>B441_face[countB441]>>B441_xs[countB441]>>B441_xe[countB441]>>B441_ys[countB441]>>B441_ye[countB441]>>B441_zs[countB441]>>B441_ze[countB441];
                        ++countB441;
						 clear(c,numint);
						 break;
                case 442: control>>B442_ID[countB442]>>B442_face[countB442]>>B442_xm[countB442]>>B442_ym[countB442]>>B442_zm[countB442]>>B442_r[countB442];
                        ++countB442;
						 clear(c,numint);
						 break;
				}
				break;

			case 'C': control>>numint;
				switch(numint)
				{
				case 75: control>>C75_x[countC75]>>C75_z[countC75]>>C75_a[countC75]>>C75_s[countC75]>>C75_l[countC75]>>C75_v[countC75];
						 ++countC75;
						 clear(c,numint);
						 break;
				}
				break;

		    case 'F': control>>numint;
				switch(numint)
				{

                case 70: control>>F70_xs[countF70]>>F70_xe[countF70]>>F70_ys[countF70]>>F70_ye[countF70]>>F70_zs[countF70]>>F70_ze[countF70];
                        ++countF70;
						 clear(c,numint);
						 break;
				case 71: control>>F71_xs[countF71]>>F71_xe[countF71]>>F71_ys[countF71]>>F71_ye[countF71]>>F71_zs[countF71]>>F71_ze[countF71];
                        ++countF71;
						 clear(c,numint);
						 break;
				case 72: control>>F72_xs[countF72]>>F72_xe[countF72]>>F72_ys[countF72]>>F72_ye[countF72]>>F72_h[countF72];
                        ++countF72;
						 clear(c,numint);
						 break;
                case 369: control>>F369_x[countF369]>>F369_z[countF369]>>F369_a[countF369]>>F369_s[countF369]>>F369_l[countF369]>>F369_v[countF369];
                        ++countF369;
						 clear(c,numint);
						 break;
				case 370: control>>F370_xs[countF370]>>F370_xe[countF370]>>F370_ys[countF370]>>F370_ye[countF370]>>F370_zs[countF370]>>F370_ze[countF370];
                        ++countF370;
						 clear(c,numint);
						 break;
				case 371: control>>F371_xs[countF371]>>F371_xe[countF371]>>F371_ys[countF371]>>F371_ye[countF371]>>F371_zs[countF371]>>F371_ze[countF371];
                        ++countF371;
						 clear(c,numint);
						 break;
               case 374: control>>F374_xc[countF374]>>F374_zc[countF374]>>F374_r[countF374];
                        ++countF374;
						 clear(c,numint);
						 break;
               case 375: control>>F375_xc[countF375]>>F375_zc[countF375]>>F375_r[countF375];
                        ++countF375;
						 clear(c,numint);
						 break;
               case 378: control>>F378_xc[countF378]>>F378_yc[countF378]>>F378_zc[countF378]>>F378_r[countF378];
                        ++countF378;
						 clear(c,numint);
						 break;
               case 379: control>>F379_xc[countF379]>>F379_yc[countF379]>>F379_zc[countF379]>>F379_r[countF379];
                        ++countF379;
						 clear(c,numint);
						 break;
				case 390: control>>F390_xs[countF390]>>F390_xe[countF390]>>F390_ys[countF390]>>F390_ye[countF390]>>F390_zs[countF390]>>F390_ze[countF390];
                        ++countF390;
						 clear(c,numint);
						 break;
				case 391: control>>F391_xs[countF391]>>F391_xe[countF391]>>F391_ys[countF391]>>F391_ye[countF391]>>F391_zs[countF391]>>F391_ze[countF391];
                        ++countF391;
						 clear(c,numint);
						 break;
               case 394: control>>F394_xc[countF394]>>F394_zc[countF394]>>F394_r[countF394];
                        ++countF394;
						 clear(c,numint);
						 break;
               case 395: control>>F395_xc[countF395]>>F395_zc[countF395]>>F395_r[countF395];
                        ++countF395;
						 clear(c,numint);
						 break;
               case 398: control>>F398_xc[countF398]>>F398_yc[countF398]>>F398_zc[countF398]>>F398_r[countF398];
                        ++countF398;
						 clear(c,numint);
						 break;
               case 399: control>>F399_xc[countF399]>>F399_yc[countF399]>>F399_zc[countF399]>>F399_r[countF399];
                        ++countF399;
						 clear(c,numint);
						 break;
				}
				break;



		    case 'P': control>>numint;
				switch(numint)
				{
				 case 35: control>>P35_ts[countP35]>>P35_te[countP35]>>P35_dt[countP35];
                        ++countP35;
						 clear(c,numint);
						 break;
                case 50: control>>P50_x[countP50]>>P50_y[countP50];
                        ++countP50;
						 clear(c,numint);
						 break;
				case 51: control>>P51_x[countP51]>>P51_y[countP51];
                        ++countP51;
						 clear(c,numint);
						 break;
                case 52: control>>P52_y[countP52];
                        ++countP52;
						 clear(c,numint);
						 break;
				case 56: control>>P56_x[countP56];
                        ++countP56;
						 clear(c,numint);
						 break;
                 case 58: control>>P58_x[countP58]>>P58_y[countP58]>>P58_T[countP58];
                        ++countP58;
						 clear(c,numint);
						 break;
				case 61: control>>P61_x[countP61]>>P61_y[countP61]>>P61_z[countP61];
                        ++countP61;
						 clear(c,numint);
						 break;
				case 62: control>>P62_xs[countP62]>>P62_xe[countP62]>>P62_ys[countP62]>>P62_ye[countP62]>>P62_zs[countP62]>>P62_ze[countP62];
                        ++countP62;
						 clear(c,numint);
						 break;
               case 63: control>>P63_x[countP63]>>P63_y[countP63];
                        ++countP63;
						 clear(c,numint);
						 break;
                case 64: control>>P64_x[countP64]>>P64_y[countP64]>>P64_z[countP64];
                        ++countP64;
						 clear(c,numint);
						 break;
				case 67: control>>P67_x[countP67];
                        ++countP67;
						 clear(c,numint);
						 break;
                case 68: control>>P68_x[countP68]>>P68_zs[countP68]>>P68_ze[countP68];
                        ++countP68;
						 clear(c,numint);
						 break;
				case 81: control>>P81_xs[countP81]>>P81_xe[countP81]>>P81_ys[countP81]>>P81_ye[countP81]>>P81_zs[countP81]>>P81_ze[countP81];
                        ++countP81;
						 clear(c,numint);
						 break;
				case 85: control>>P85_x[countP85]>>P85_y[countP85]>>P85_r[countP85]>>P85_cd[countP85]>>P85_cm[countP85];
                        ++countP85;
						 clear(c,numint);
						 break;
				case 121: control>>P121_x[countP121]>>P121_y[countP121];
                        ++countP121;
						 clear(c,numint);
						 break;
				case 123: control>>P123_y[countP123];
                        ++countP123;
						 clear(c,numint);
						 break;
				case 124: control>>P124_x[countP124];
                        ++countP124;
						 clear(c,numint);
						 break;
				case 125: control>>P125_x[countP125]>>P125_y[countP125];
                        ++countP125;
						 clear(c,numint);
						 break;
               case 184: control>>P184_its[countP184]>>P184_ite[countP184]>>P184_dit[countP184];
                        ++countP184;
						 clear(c,numint);
						 break;
               case 185: control>>P185_ts[countP185]>>P185_te[countP185]>>P185_dt[countP185];
                        ++countP185;
						 clear(c,numint);
						 break;
               case 194: control>>P194_its[countP194]>>P194_ite[countP194]>>P194_dit[countP194];
                        ++countP194;
						 clear(c,numint);
						 break;
               case 195: control>>P195_ts[countP195]>>P195_te[countP195]>>P195_dt[countP195];
                        ++countP195;
						 clear(c,numint);
						 break;
               case 230: control>>P230_x[countP230];
                        ++countP230;
						 clear(c,numint);
						 break;
               case 240: control>>P240_x[countP240];
                        ++countP240;
						 clear(c,numint);
						 break;
				case 351: control>>P351_x[countP351]>>P351_y[countP351];
                        ++countP351;
						 clear(c,numint);
						 break;
				case 352: control>>P352_x[countP352]>>P352_y[countP352];
                        ++countP352;
						 clear(c,numint);
						 break;
				}
				break;
                
            case 'Q': control>>numint;
				switch(numint)
				{

                case 110: control>>Q110_xs[countQ110]>>Q110_xe[countQ110]>>Q110_ys[countQ110]>>Q110_ye[countQ110]>>Q110_zs[countQ110]>>Q110_ze[countQ110];
                        ++countQ110;
						 clear(c,numint);
						 break;
				}
				break;

            case 'S': control>>numint;
				switch(numint)
				{
				case 73: control>>S73_val[countS73]>>S73_dist[countS73]>>S73_b[countS73]>>S73_x[countS73]>>S73_y[countS73];
                        ++countS73;
						 clear(c,numint);
						 break;
				}
				break;

            case 'W': control>>numint;
				switch(numint)
				{
				case 41: control>>W41_xc[countW41]>>W41_yc[countW41]>>W41_zs[countW41]>>W41_ze[countW41]>>W41_vel[countW41]>>W41_beta[countW41];
                        ++countW41;
						 clear(c,numint);
						 break;
				}
				break;


			case 'X': control>>numint;
				switch(numint)
				{

                case 110: control>>X110_xs[countX110]>>X110_xe[countX110]>>X110_ys[countX110]>>X110_ye[countX110]>>X110_zs[countX110]>>X110_ze[countX110];
                        ++countX110;
						 clear(c,numint);
						 break;
                case 163: control>>X163_x1[countX163]>>X163_y1[countX163]>>X163_z1[countX163]>>X163_x2[countX163]>>X163_y2[countX163]>>X163_z2[countX163]
                                 >>X163_x3[countX163]>>X163_y3[countX163]>>X163_z3[countX163]>>X163_x4[countX163]>>X163_y4[countX163]>>X163_z4[countX163]
                                 >>X163_x5[countX163]>>X163_y5[countX163]>>X163_z5[countX163]>>X163_x6[countX163]>>X163_y6[countX163]>>X163_z6[countX163];
                        ++countX163;
						 clear(c,numint);
						 break;
                case 164: control>>X164_x1[countX164]>>X164_y1[countX164]>>X164_z1[countX164]>>X164_x2[countX164]>>X164_y2[countX164]>>X164_z2[countX164]
                                 >>X164_x3[countX164]>>X164_y3[countX164]>>X164_z3[countX164]>>X164_x4[countX164]>>X164_y4[countX164]>>X164_z4[countX164]
                                 >>X164_x5[countX164]>>X164_y5[countX164]>>X164_z5[countX164]>>X164_x6[countX164]>>X164_y6[countX164]>>X164_z6[countX164]
                                 >>X164_x7[countX164]>>X164_y7[countX164]>>X164_z7[countX164]>>X164_x8[countX164]>>X164_y8[countX164]>>X164_z8[countX164];
                        ++countX164;
						 clear(c,numint);
						 break;
				case 311: control>>X311_xs[countX311]>>X311_xe[countX311]>>X311_ys[countX311]>>X311_ye[countX311]>>X311_zs[countX311]>>X311_ze[countX311]
								>>X311_w[countX311]>>X311_rho_c[countX311]>>X311_EA[countX311]>>X311_d[countX311]>>X311_l[countX311]>>X311_H[countX311]
								>>X311_P[countX311]>>X311_facT[countX311];
                        ++countX311;
						 clear(c,numint);
						 break;
                case 312: control>>X311_xs[countX312]>>X311_xe[countX312]>>X311_ys[countX312]>>X311_ye[countX312]>>X311_zs[countX312]>>X311_ze[countX312]
								>>X312_k[countX312]>>X312_T0[countX312];
                        ++countX312;
						 clear(c,numint);
						 break;
                case 314: 
                         for (int i = 0; i < mooring_count; i++)
                              control>>X314_T[i];
						 
                         clear(c,numint);
						 break;
                case 315: 
                         for (int i = 0; i < mooring_count; i++)
                              control>>X315_t[i];
						 
                         clear(c,numint);
						 break;
				case 320: control>>X320_type[countX320];
                        ++countX320;
						 clear(c,numint);
						 break;
				case 321: control>>X321_Sn[countX321]>>X321_d[countX321]>>X321_lambda[countX321]>>X321_dk[countX321]>>X321_rho[countX321]>>X321_nd[countX321]>>X321_nl[countX321];
                        ++countX321;
						 clear(c,numint);
						 break;
				case 322: control>>X322_D[countX322]>>X322_L[countX322]>>X322_x0[countX322]>>X322_y0[countX322]>>X322_z0[countX322]>>X322_phi[countX322]>>X322_theta[countX322]>>X322_psi[countX322];
                        ++countX322;
						 clear(c,numint);
						 break;
				case 324: control>>X324_x[countX324]>>X324_y[countX324]>>X324_z[countX324];
                        ++countX324;
						 clear(c,numint);
						 break;
				}
				break;
				
			case 'Z': control>>numint;
				switch(numint)
				{
                case 11: control>>Z11_x[countZ11]>>Z11_y[countZ11]>>Z11_z[countZ11]>>Z11_l[countZ11]>>Z11_w[countZ11]>>Z11_t[countZ11]>>Z11_rho[countZ11]>>Z11_e[countZ11]>>Z11_ix[countZ11]>>Z11_iy[countZ11]>>Z11_iz[countZ11]>>Z11_nu[countZ11]>>Z11_n[countZ11];
                        ++countZ11;
						 clear(c,numint);
						 break;
                case 12: control>>Z12_cdx>>Z12_cdy>>Z12_cdz>>Z12_ckx>>Z12_cky>>Z12_ckz;
						 clear(c,numint);
						 break;
				}
				break;
			}
		}
        if(count>1e7)
        {
        cout<<endl;
        cout<<"!!! missing input parameter in ctrl.txt !!!"<<endl<<endl;
        cout<<"!!! please check the REEF3D User Guide !!!"<<endl<<endl<<endl<<endl;

        exit(0);
        }
	}

	control.close();
}

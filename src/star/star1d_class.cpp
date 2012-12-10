#include"star.h"
#include<string.h>
#include<stdlib.h>

star1d::star1d() {}

star1d::~star1d() {
}

star1d::star1d(const star1d &A) : star2d(A) {}

star1d &star1d::operator=(const star1d &A) {

	star2d::operator=(A);
	
	return *this;

}

void star1d::write_tag(OUTFILE *fp,char mode) const {

	char tag[7]="star1d";
	
	if(mode=='b') fp->write("tag",tag,7);
	else fp->write_fmt("tag","%s",&tag);
		
}

int star1d::check_tag(const char *tag) const {

	if(strcmp(tag,"star1d")) return 0;
	return 1;

}

int star1d::read_old(const char *input_file){

	FILE *fp;
	char tag[7],mode,*c;
	int ndom,i;
	
	if(!(fp=fopen(input_file,"rb"))) {
		return 0;
	}
	fread(tag,1,7,fp);
	tag[6]='\0';
	if(strcmp(tag,"star1d")) {
		return 0;
	}
	fread(&mode,1,1,fp);
	fclose(fp);

	if(mode=='b')
		fp=fopen(input_file,"rb");
	else
		fp=fopen(input_file,"rt");
	
	fread(tag,1,7,fp);
	fread(&mode,1,1,fp);
	
	if(mode=='b') {
		fread(&ndom,sizeof(int),1,fp);
		map.gl.set_ndomains(ndom);
		fread(map.gl.npts,sizeof(int),ndom,fp);
		fread(map.gl.xif,sizeof(double),ndom+1,fp);
		fread(&M,sizeof(double),1,fp);
		fread(&R,sizeof(double),1,fp);
		fread(&X,sizeof(double),1,fp);
		fread(&Z,sizeof(double),1,fp);
		fread(&Xc,sizeof(double),1,fp);
		fread(&conv,sizeof(int),1,fp);
		fread(&surff,sizeof(double),1,fp);
		fread(&Tc,sizeof(double),1,fp);
		fread(&pc,sizeof(double),1,fp);
		c=opa.name;
		*c=fgetc(fp);
		while(*(c++)!='\0') *c=fgetc(fp);
		c=eos.name;
		*c=fgetc(fp);
		while(*(c++)!='\0') *c=fgetc(fp);
		c=nuc.name;
		*c=fgetc(fp);
		while(*(c++)!='\0') *c=fgetc(fp);
		c=atm_name;
		*c=fgetc(fp);
		while(*(c++)!='\0') *c=fgetc(fp);
	} else {
		fscanf(fp,"\n%d ",&ndom);
		map.gl.set_ndomains(ndom);
		for(i=0;i<ndom;i++) fscanf(fp,"%d ",(map.gl.npts+i));
		for(i=0;i<ndom+1;i++) fscanf(fp,"%le ",(map.gl.xif+i));
		fscanf(fp,"\n%le %le %le %le\n",&M,&R,&X,&Z);
		fscanf(fp,"%le %d %le\n",&Xc,&conv,&surff);		
		fscanf(fp,"%le %le\n",&Tc,&pc);
		fscanf(fp,"%s\n",opa.name);
		fscanf(fp,"%s\n",eos.name);
		fscanf(fp,"%s\n",nuc.name);
		fscanf(fp,"%s\n",atm_name);
	}
	map.leg.npts=1;
	map.init();
	
	phi.read(nr,1,fp,mode);
	p.read(nr,1,fp,mode);
	T.read(nr,1,fp,mode);
	
	fclose(fp);
	fill();
	return 1;
		
}

int star1d::init(const char *input_file,const char *param_file,int argc,char *argv[]) {

	cmdline_parser cmd;
	file_parser fp;
	char *arg,*val,default_params[256];
	mapping map0;
	int i,k,change_grid=0;
	matrix Tr;

	sprintf(default_params,"%s/config/1d_default.par",ESTER_ROOT);

	if(*input_file) {
		if(!read(input_file)) {
			printf("Error reading input file: %s\n",input_file);
			return 0;
		}
		map0=map;
	} else {
		if(!fp.open(default_params)) { 
			printf("Can't open default parameters file %s\n",default_params);
			return 0;
		}
		else {
			while(k=fp.get(arg,val)) {			
				if(i=check_arg(arg,val,&change_grid)) {
					printf("Sintax error in parameters file %s, line %d\n",param_file,k);
					if(i==2) {
						printf("Error: Argument to '%s' missing\n",arg);
						return 0;
					}
					if(i==1) {
						printf("Unknown parameter %s\n",arg);
						return 0;
					}
				}
			}
			fp.close();
		}
		change_grid=0;
	}
	
	if(*param_file) {
		if(!fp.open(param_file)) { 
			printf("Can't open parameters file %s\n",param_file);
			return 0;
		}
		else {
			while(k=fp.get(arg,val)) {			
				if(i=check_arg(arg,val,&change_grid)) {
					printf("Sintax error in parameters file %s, line %d\n",param_file,k);
					if(i==2) {
						printf("Error: Argument to '%s' missing\n",arg);
						return 0;
					}
					if(i==1) {
						printf("Unknown parameter %s\n",arg);
						return 0;
					}
				}
			}
			fp.close();
		}
	}
	
	cmd.open(argc,argv);
	while(int err_code=cmd.get(arg,val)) {
		if(err_code==-1) exit(1);
		err_code=check_arg(arg,val,&change_grid);
		if(err_code==2) {
			fprintf(stderr,"Error: Argument to '%s' missing\n",arg);
			return 0;
		}
		if(err_code==1) {
			fprintf(stderr,"Unknown parameter '%s'\n",arg);
			return 0;
		}
		cmd.ack(arg,val);
	}
	cmd.close();
	
	if((change_grid&1)&&!(change_grid&2)) {
		fprintf(stderr,"Must specify number of points per domain (npts)\n");
		exit(1);
	}
	if (*input_file) {
		if(change_grid) {
			mapping map_new;
			map_new=map;
			map=map0;
			remap(map_new.ndomains,map_new.gl.npts,map_new.nth,map_new.nex);
		}
	} else {
		map.leg.npts=1;
		for(i=0;i<=ndomains;i++) map.gl.xif[i]=i*1./ndomains;
		map.init();
		T=1-0.5*r*r;
		p=T;
		phi=-T;
		G=0*T;
		w=0*T;
		conv=0;
		phiex=zeros(map.nex,map.nth);
	}
	
	fill();
	
	return 1;
	
}

int star1d::check_arg(char *arg,char *val,int *change_grid) {

	int err=0,i;
	char *tok;

	if(!strcmp(arg,"ndomains")) {
		if(val==NULL) return 2;
		map.gl.set_ndomains(atoi(val));
		*change_grid=*change_grid|1;
		if(*change_grid&2) {
			fprintf(stderr,"ndomains must be specified before npts\n");
			exit(1);
		}
	}
	else if(!strcmp(arg,"npts")) {
		if(val==NULL) return 2;
		tok=strtok(val,",");
		i=0;
		while(tok!=NULL) {
			*(map.gl.npts+i)=atoi(tok);
			tok=strtok(NULL,",");
			i++;
		}
		if(i==1) {
			for(i=1;i<map.gl.ndomains;i++) {
				*(map.gl.npts+i)=*map.gl.npts;
			}
		}	
		*change_grid=*change_grid|2;
	}
	else if(!strcmp(arg,"xif")) {
		fprintf(stderr,"Warning: Parameter xif is now automatically handled by the code and cannot be modified by the user\n"); 
	}
	else if(!strcmp(arg,"M")) {
		if(val==NULL) return 2;
		M=atof(val)*M_SUN;
	}
	else if(!strcmp(arg,"X")) {
		if(val==NULL) return 2;
		X=atof(val);
	}
	else if(!strcmp(arg,"Z")) {
		if(val==NULL) return 2;
		Z=atof(val);
	}
	else if(!strcmp(arg,"Xc")) {
		if(val==NULL) return 2;
		Xc=atof(val);
	}
	else if(!strcmp(arg,"conv")) {
		fprintf(stderr,"Param. conv is no longer modifiable. Disable core convection with core_convec 0.\n");
	}
	else if(!strcmp(arg,"surff")) {
		if(val==NULL) return 2;
		surff=atof(val);
	}
	else if(!strcmp(arg,"Tc")) {
		if(val==NULL) return 2;
		Tc=atof(val);
	}
	else if(!strcmp(arg,"pc")) {
		if(val==NULL) return 2;
		pc=atof(val);
	}
	else if(!strcmp(arg,"opa")) {
		if(val==NULL) return 2;
		strcpy(opa.name,val);
	}
	else if(!strcmp(arg,"eos")) {
		if(val==NULL) return 2;
		strcpy(eos.name,val);
	}
	else if(!strcmp(arg,"nuc")) {
		if(val==NULL) return 2;
		strcpy(nuc.name,val);
	}
	else if(!strcmp(arg,"atm")) {
		if(val==NULL) return 2;
		strcpy(atm_name,val);
	}
	else err=1;

	return err;

}



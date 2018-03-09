#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define _set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define _get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)

int main(int argc, char** argv) {

  if (argc < 3) {
    printf("USAGE: %s refgenome fpga_pac\n", argv[0]);
    return 0;
  }

  char pac_path[200];
  sprintf(pac_path, "%s.pac", argv[1]);
  FILE* pac_file = fopen(pac_path, "rb");
  FILE* pac_fpga = fopen(argv[2], "rb");

  if (!pac_file) {
    printf("cannot find pac file: %s\n", pac_path);
    return 1;
  }

  if (!pac_fpga) {
    printf("cannot find fpga.pac file: %s\n", pac_fpga);
    return 1;
  }

  int pac_size_1 = 0;
  int pac_size_2 = 0;

  fseek(pac_file, 0, SEEK_END);
  pac_size_1 = ftell(pac_file);
  fseek(pac_file, 0, SEEK_SET);

  fseek(pac_fpga, 0, SEEK_END);
  pac_size_2 = ftell(pac_fpga);
  fseek(pac_fpga, 0, SEEK_SET);

  if (pac_size_1*2-1 != pac_size_2) {
    printf("size mismatch\n");
    return 2;
  }


  char* pac_1 = (char*)calloc(pac_size_2, sizeof(char));
  char* pac_2 = (char*)calloc(pac_size_2, sizeof(char));

  int remain = 0;

  fread(pac_1, 1, pac_size_1-1, pac_file);
  fread(&remain, 1, 1, pac_file);
  fclose(pac_file);

  fread(pac_2, 1, pac_size_2, pac_fpga);
  fclose(pac_fpga);

  // the last element stores the remain

  int64_t l_pac = ((int64_t)pac_size_1-2)*4 + remain;

  printf("size 1 = %ld\n", pac_size_1);
  printf("size 2 = %ld\n", pac_size_2);
  printf("l_pac = %ld\n", l_pac);
  printf("remain = %d\n", remain);

  int64_t k = l_pac;
  for (int64_t l = l_pac - 1; l >= 0; --l) {
    _set_pac(pac_1, k, 3-_get_pac(pac_1, l));
    k++;
  }

  for (int k = 0; k < pac_size_2-1; k++) {
    if (pac_1[k] != pac_2[k]) {
      printf("element mismatch at %d: %d != %d\n", k, pac_1[k], pac_2[k]);
      break;
    }
  }

  printf("matched\n");

  free(pac_1);
  free(pac_2);

  return 0;
}

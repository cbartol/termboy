#include <gb/gb.hpp>
#include <fstream>
#include <ctime>

#define CPU_CPP
namespace GameBoy {


#include "mmio.cpp"
#include "memory.cpp"
#include "timing.cpp"
#include "serialization.cpp"
CPU cpu;

void CPU::op_xx() {
}

void CPU::op_cb() {
  exec_cb();
}

//8-bit load commands

template<unsigned x, unsigned y> void CPU::op_ld_r_r() {
  r[x] = r[y];
}

template<unsigned x> void CPU::op_ld_r_n() {
  r[x] = op_read(r[PC]++);
}

template<unsigned x> void CPU::op_ld_r_hl() {
  r[x] = op_read(r[HL]);
}

template<unsigned x> void CPU::op_ld_hl_r() {
  op_write(r[HL], r[x]);
}

void CPU::op_ld_hl_n() {
  op_write(r[HL], op_read(r[PC]++));
}

template<unsigned x> void CPU::op_ld_a_rr() {
  r[A] = op_read(r[x]);
}

void CPU::op_ld_a_nn() {
  uint8 lo = op_read(r[PC]++);
  uint8 hi = op_read(r[PC]++);
  r[A] = op_read((hi << 8) | (lo << 0));
}

template<unsigned x> void CPU::op_ld_rr_a() {
  op_write(r[x], r[A]);
}

void CPU::op_ld_nn_a() {
  uint8 lo = op_read(r[PC]++);
  uint8 hi = op_read(r[PC]++);
  op_write((hi << 8) | (lo << 0), r[A]);
}

void CPU::op_ld_a_ffn() {
  r[A] = op_read(0xff00 + op_read(r[PC]++));
}

void CPU::op_ld_ffn_a() {
  op_write(0xff00 + op_read(r[PC]++), r[A]);
}

void CPU::op_ld_a_ffc() {
  r[A] = op_read(0xff00 + r[C]);
}

void CPU::op_ld_ffc_a() {
  op_write(0xff00 + r[C], r[A]);
}

void CPU::op_ldi_hl_a() {
  op_write(r[HL], r[A]);
  r[HL]++;
}

void CPU::op_ldi_a_hl() {
  r[A] = op_read(r[HL]);
  r[HL]++;
}

void CPU::op_ldd_hl_a() {
  op_write(r[HL], r[A]);
  r[HL]--;
}

void CPU::op_ldd_a_hl() {
  r[A] = op_read(r[HL]);
  r[HL]--;
}

//16-bit load commands

template<unsigned x> void CPU::op_ld_rr_nn() {
  r[x]  = op_read(r[PC]++) << 0;
  r[x] |= op_read(r[PC]++) << 8;
}

void CPU::op_ld_nn_sp() {
  uint16 addr = op_read(r[PC]++) << 0;
  addr |= op_read(r[PC]++) << 8;
  op_write(addr + 0, r[SP] >> 0);
  op_write(addr + 1, r[SP] >> 8);
}

void CPU::op_ld_sp_hl() {
  r[SP] = r[HL];
  op_io();
}

template<unsigned x> void CPU::op_push_rr() {
  op_write(--r[SP], r[x] >> 8);
  op_write(--r[SP], r[x] >> 0);
  op_io();
}

template<unsigned x> void CPU::op_pop_rr() {
  r[x]  = op_read(r[SP]++) << 0;
  r[x] |= op_read(r[SP]++) << 8;
}

//8-bit arithmetic commands

void CPU::opi_add_a(uint8 x) {
  uint16 rh = r[A] + x;
  uint16 rl = (r[A] & 0x0f) + (x & 0x0f);
  r[A] = rh;
  r.f.z = (uint8)rh == 0;
  r.f.n = 0;
  r.f.h = rl > 0x0f;
  r.f.c = rh > 0xff;
}

template<unsigned x> void CPU::op_add_a_r() { opi_add_a(r[x]); }
void CPU::op_add_a_n() { opi_add_a(op_read(r[PC]++)); }
void CPU::op_add_a_hl() { opi_add_a(op_read(r[HL])); }

void CPU::opi_adc_a(uint8 x) {
  uint16 rh = r[A] + x + r.f.c;
  uint16 rl = (r[A] & 0x0f) + (x & 0x0f) + r.f.c;
  r[A] = rh;
  r.f.z = (uint8)rh == 0;
  r.f.n = 0;
  r.f.h = rl > 0x0f;
  r.f.c = rh > 0xff;
}

template<unsigned x> void CPU::op_adc_a_r() { opi_adc_a(r[x]); }
void CPU::op_adc_a_n() { opi_adc_a(op_read(r[PC]++)); }
void CPU::op_adc_a_hl() { opi_adc_a(op_read(r[HL])); }

void CPU::opi_sub_a(uint8 x) {
  uint16 rh = r[A] - x;
  uint16 rl = (r[A] & 0x0f) - (x & 0x0f);
  r[A] = rh;
  r.f.z = (uint8)rh == 0;
  r.f.n = 1;
  r.f.h = rl > 0x0f;
  r.f.c = rh > 0xff;
}

template<unsigned x> void CPU::op_sub_a_r() { opi_sub_a(r[x]); }
void CPU::op_sub_a_n() { opi_sub_a(op_read(r[PC]++)); }
void CPU::op_sub_a_hl() { opi_sub_a(op_read(r[HL])); }

void CPU::opi_sbc_a(uint8 x) {
  uint16 rh = r[A] - x - r.f.c;
  uint16 rl = (r[A] & 0x0f) - (x & 0x0f) - r.f.c;
  r[A] = rh;
  r.f.z = (uint8)rh == 0;
  r.f.n = 1;
  r.f.h = rl > 0x0f;
  r.f.c = rh > 0xff;
}

template<unsigned x> void CPU::op_sbc_a_r() { opi_sbc_a(r[x]); }
void CPU::op_sbc_a_n() { opi_sbc_a(op_read(r[PC]++)); }
void CPU::op_sbc_a_hl() { opi_sbc_a(op_read(r[HL])); }

void CPU::opi_and_a(uint8 x) {
  r[A] &= x;
  r.f.z = r[A] == 0;
  r.f.n = 0;
  r.f.h = 1;
  r.f.c = 0;
}

template<unsigned x> void CPU::op_and_a_r() { opi_and_a(r[x]); }
void CPU::op_and_a_n() { opi_and_a(op_read(r[PC]++)); }
void CPU::op_and_a_hl() { opi_and_a(op_read(r[HL])); }

void CPU::opi_xor_a(uint8 x) {
  r[A] ^= x;
  r.f.z = r[A] == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = 0;
}

template<unsigned x> void CPU::op_xor_a_r() { opi_xor_a(r[x]); }
void CPU::op_xor_a_n() { opi_xor_a(op_read(r[PC]++)); }
void CPU::op_xor_a_hl() { opi_xor_a(op_read(r[HL])); }

void CPU::opi_or_a(uint8 x) {
  r[A] |= x;
  r.f.z = r[A] == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = 0;
}

template<unsigned x> void CPU::op_or_a_r() { opi_or_a(r[x]); }
void CPU::op_or_a_n() { opi_or_a(op_read(r[PC]++)); }
void CPU::op_or_a_hl() { opi_or_a(op_read(r[HL])); }

void CPU::opi_cp_a(uint8 x) {
  uint16 rh = r[A] - x;
  uint16 rl = (r[A] & 0x0f) - (x & 0x0f);
  r.f.z = (uint8)rh == 0;
  r.f.n = 1;
  r.f.h = rl > 0x0f;
  r.f.c = rh > 0xff;
}

template<unsigned x> void CPU::op_cp_a_r() { opi_cp_a(r[x]); }
void CPU::op_cp_a_n() { opi_cp_a(op_read(r[PC]++)); }
void CPU::op_cp_a_hl() { opi_cp_a(op_read(r[HL])); }

template<unsigned x> void CPU::op_inc_r() {
  r[x]++;
  r.f.z = r[x] == 0;
  r.f.n = 0;
  r.f.h = (r[x] & 0x0f) == 0x00;
}

void CPU::op_inc_hl() {
  uint8 n = op_read(r[HL]);
  op_write(r[HL], ++n);
  r.f.z = n == 0;
  r.f.n = 0;
  r.f.h = (n & 0x0f) == 0x00;
}

template<unsigned x> void CPU::op_dec_r() {
  r[x]--;
  r.f.z = r[x] == 0;
  r.f.n = 1;
  r.f.h = (r[x] & 0x0f) == 0x0f;
}

void CPU::op_dec_hl() {
  uint8 n = op_read(r[HL]);
  op_write(r[HL], --n);
  r.f.z = n == 0;
  r.f.n = 1;
  r.f.h = (n & 0x0f) == 0x0f;
}

void CPU::op_daa() {
  uint16 a = r[A];
  if(r.f.n == 0) {
    if(r.f.h || (a & 0x0f) > 0x09) a += 0x06;
    if(r.f.c || (a       ) > 0x9f) a += 0x60;
  } else {
    if(r.f.h) {
      a -= 0x06;
      if(r.f.c == 0) a &= 0xff;
    }
    if(r.f.c) a -= 0x60;
  }
  r[A] = a;
  r.f.z = r[A] == 0;
  r.f.h = 0;
  r.f.c |= a & 0x100;
}

void CPU::op_cpl() {
  r[A] ^= 0xff;
  r.f.n = 1;
  r.f.h = 1;
}

//16-bit arithmetic commands

template<unsigned x> void CPU::op_add_hl_rr() {
  op_io();
  uint32 rb = (r[HL] + r[x]);
  uint32 rn = (r[HL] & 0xfff) + (r[x] & 0xfff);
  r[HL] = rb;
  r.f.n = 0;
  r.f.h = rn > 0x0fff;
  r.f.c = rb > 0xffff;
}

template<unsigned x> void CPU::op_inc_rr() {
  op_io();
  r[x]++;
}

template<unsigned x> void CPU::op_dec_rr() {
  op_io();
  r[x]--;
}

void CPU::op_add_sp_n() {
  op_io();
  op_io();
  signed n = (int8)op_read(r[PC]++);
  r.f.z = 0;
  r.f.n = 0;
  r.f.h = ((r[SP] & 0x0f) + (n & 0x0f)) > 0x0f;
  r.f.c = ((r[SP] & 0xff) + (n & 0xff)) > 0xff;
  r[SP] += n;
}

void CPU::op_ld_hl_sp_n() {
  op_io();
  signed n = (int8)op_read(r[PC]++);
  r.f.z = 0;
  r.f.n = 0;
  r.f.h = ((r[SP] & 0x0f) + (n & 0x0f)) > 0x0f;
  r.f.c = ((r[SP] & 0xff) + (n & 0xff)) > 0xff;
  r[HL] = r[SP] + n;
}

//rotate/shift commands

void CPU::op_rlca() {
  r[A] = (r[A] << 1) | (r[A] >> 7);
  r.f.z = 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = r[A] & 0x01;
}

void CPU::op_rla() {
  bool c = r[A] & 0x80;
  r[A] = (r[A] << 1) | (r.f.c << 0);
  r.f.z = 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = c;
}

void CPU::op_rrca() {
  r[A] = (r[A] >> 1) | (r[A] << 7);
  r.f.z = 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = r[A] & 0x80;
}

void CPU::op_rra() {
  bool c = r[A] & 0x01;
  r[A] = (r[A] >> 1) | (r.f.c << 7);
  r.f.z = 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = c;
}

template<unsigned x> void CPU::op_rlc_r() {
  r[x] = (r[x] << 1) | (r[x] >> 7);
  r.f.z = r[x] == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = r[x] & 0x01;
}

void CPU::op_rlc_hl() {
  uint8 n = op_read(r[HL]);
  n = (n << 1) | (n >> 7);
  op_write(r[HL], n);
  r.f.z = n == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = n & 0x01;
}

template<unsigned x> void CPU::op_rl_r() {
  bool c = r[x] & 0x80;
  r[x] = (r[x] << 1) | (r.f.c << 0);
  r.f.z = r[x] == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = c;
}

void CPU::op_rl_hl() {
  uint8 n = op_read(r[HL]);
  bool c = n & 0x80;
  n = (n << 1) | (r.f.c << 0);
  op_write(r[HL], n);
  r.f.z = n == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = c;
}

template<unsigned x> void CPU::op_rrc_r() {
  r[x] = (r[x] >> 1) | (r[x] << 7);
  r.f.z = r[x] == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = r[x] & 0x80;
}

void CPU::op_rrc_hl() {
  uint8 n = op_read(r[HL]);
  n = (n >> 1) | (n << 7);
  op_write(r[HL], n);
  r.f.z = n == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = n & 0x80;
}

template<unsigned x> void CPU::op_rr_r() {
  bool c = r[x] & 0x01;
  r[x] = (r[x] >> 1) | (r.f.c << 7);
  r.f.z = r[x] == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = c;
}

void CPU::op_rr_hl() {
  uint8 n = op_read(r[HL]);
  bool c = n & 0x01;
  n = (n >> 1) | (r.f.c << 7);
  op_write(r[HL], n);
  r.f.z = n == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = c;
}

template<unsigned x> void CPU::op_sla_r() {
  bool c = r[x] & 0x80;
  r[x] <<= 1;
  r.f.z = r[x] == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = c;
}

void CPU::op_sla_hl() {
  uint8 n = op_read(r[HL]);
  bool c = n & 0x80;
  n <<= 1;
  op_write(r[HL], n);
  r.f.z = n == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = c;
}

template<unsigned x> void CPU::op_swap_r() {
  r[x] = (r[x] << 4) | (r[x] >> 4);
  r.f.z = r[x] == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = 0;
}

void CPU::op_swap_hl() {
  uint8 n = op_read(r[HL]);
  n = (n << 4) | (n >> 4);
  op_write(r[HL], n);
  r.f.z = n == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = 0;
}

template<unsigned x> void CPU::op_sra_r() {
  bool c = r[x] & 0x01;
  r[x] = (int8)r[x] >> 1;
  r.f.z = r[x] == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = c;
}

void CPU::op_sra_hl() {
  uint8 n = op_read(r[HL]);
  bool c = n & 0x01;
  n = (int8)n >> 1;
  op_write(r[HL], n);
  r.f.z = n == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = c;
}

template<unsigned x> void CPU::op_srl_r() {
  bool c = r[x] & 0x01;
  r[x] >>= 1;
  r.f.z = r[x] == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = c;
}

void CPU::op_srl_hl() {
  uint8 n = op_read(r[HL]);
  bool c = n & 0x01;
  n >>= 1;
  op_write(r[HL], n);
  r.f.z = n == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = c;
}

//single-bit commands

template<unsigned b, unsigned x> void CPU::op_bit_n_r() {
  r.f.z = (r[x] & (1 << b)) == 0;
  r.f.n = 0;
  r.f.h = 1;
}

template<unsigned b> void CPU::op_bit_n_hl() {
  uint8 n = op_read(r[HL]);
  r.f.z = (n & (1 << b)) == 0;
  r.f.n = 0;
  r.f.h = 1;
}

template<unsigned b, unsigned x> void CPU::op_set_n_r() {
  r[x] |= 1 << b;
}

template<unsigned b> void CPU::op_set_n_hl() {
  uint8 n = op_read(r[HL]);
  n |= 1 << b;
  op_write(r[HL], n);
}

template<unsigned b, unsigned x> void CPU::op_res_n_r() {
  r[x] &= ~(1 << b);
}

template<unsigned b> void CPU::op_res_n_hl() {
  uint8 n = op_read(r[HL]);
  n &= ~(1 << b);
  op_write(r[HL], n);
}

//control commands

void CPU::op_ccf() {
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = !r.f.c;
}

void CPU::op_scf() {
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = 1;
}

void CPU::op_nop() {
}

void CPU::op_halt() {
  r.halt = true;
  while(r.halt == true) op_io();
}

void CPU::op_stop() {
  if(stop()) return;
  r.stop = true;
  while(r.stop == true) op_io();
}

void CPU::op_di() {
  r.ime = 0;
}

void CPU::op_ei() {
  r.ei = true;
//r.ime = 1;
}

//jump commands

void CPU::op_jp_nn() {
  uint8 lo = op_read(r[PC]++);
  uint8 hi = op_read(r[PC]++);
  r[PC] = (hi << 8) | (lo << 0);
  op_io();
}

void CPU::op_jp_hl() {
  r[PC] = r[HL];
}

template<unsigned x, bool y> void CPU::op_jp_f_nn() {
  uint8 lo = op_read(r[PC]++);
  uint8 hi = op_read(r[PC]++);
  if(r.f[x] == y) {
    r[PC] = (hi << 8) | (lo << 0);
    op_io();
  }
}

void CPU::op_jr_n() {
  int8 n = op_read(r[PC]++);
  r[PC] += n;
  op_io();
}

template<unsigned x, bool y> void CPU::op_jr_f_n() {
  int8 n = op_read(r[PC]++);
  if(r.f[x] == y) {
    r[PC] += n;
    op_io();
  }
}

void CPU::op_call_nn() {
  uint8 lo = op_read(r[PC]++);
  uint8 hi = op_read(r[PC]++);
  op_write(--r[SP], r[PC] >> 8);
  op_write(--r[SP], r[PC] >> 0);
  r[PC] = (hi << 8) | (lo << 0);
  op_io();
}

template<unsigned x, bool y> void CPU::op_call_f_nn() {
  uint8 lo = op_read(r[PC]++);
  uint8 hi = op_read(r[PC]++);
  if(r.f[x] == y) {
    op_write(--r[SP], r[PC] >> 8);
    op_write(--r[SP], r[PC] >> 0);
    r[PC] = (hi << 8) | (lo << 0);
    op_io();
  }
}

void CPU::op_ret() {
  uint8 lo = op_read(r[SP]++);
  uint8 hi = op_read(r[SP]++);
  r[PC] = (hi << 8) | (lo << 0);
  op_io();
}

template<unsigned x, bool y> void CPU::op_ret_f() {
  op_io();
  if(r.f[x] == y) {
    uint8 lo = op_read(r[SP]++);
    uint8 hi = op_read(r[SP]++);
    r[PC] = (hi << 8) | (lo << 0);
    op_io();
  }
}

void CPU::op_reti() {
  uint8 lo = op_read(r[SP]++);
  uint8 hi = op_read(r[SP]++);
  r[PC] = (hi << 8) | (lo << 0);
  op_io();
  r.ime = 1;
}

template<unsigned n> void CPU::op_rst_n() {
  op_write(--r[SP], r[PC] >> 8);
  op_write(--r[SP], r[PC] >> 0);
  r[PC] = n;
  op_io();
}

void CPU::dump() {
    std::ofstream txt, csv;
    txt.open("/tmp/benchmarks.txt", std::ios::out);
    csv.open("/tmp/benchmarks.csv", std::ios::out);

    txt  << "Benchmarks for Termboy total instructions: " << instruction_count << std::endl;
    txt  << "-------------------------------------------" << std::endl;

    csv  << "INSTRUCTION, COUNT, AVG TIME, MAX TIME";


    for(int i = 0 ; i < 256 ; i++) {
        if(inst_counter[i] != 0) {
           txt    << "I-COUNT["<< std::hex << i << "]: " 
                  << std::dec << inst_counter[i] 
                  << "\tTIME: "<< std::dec << times[i]
                  << "\tMAX: "<< max_times[i] << std::endl;

           csv << std::hex << i << std:: dec << "," <<  times[i] << "," << max_times[i] << std::endl;
            }
    }

    txt  << "-------------------------------------------" << std::endl;
    txt  << " CB benchmarks "  << inst_counter[0xcb] << std::endl;
    txt  << "-------------------------------------------" << std::endl;
    for(int i = 0 ; i < 256 ; i++) {
        if(cb_inst_counter[i] != 0) {
           txt    << "I-COUNT["<< std::hex << i << "]: " 
                  << std::dec << cb_inst_counter[i] 
                  << "\tTIME: "<< std::dec << cb_times[i]
                  << "\tMAX:  "<< cb_max_times[i] << std::endl;

          // csv << std::hex << i << std:: dec << "," <<  times[i] << "," << max_times[i] << std::endl;
            }
    }

    txt.close();
    csv.close();
}

void CPU::get_times(double final_time) {
    if(max_times[last_inst] < final_time )
        max_times[last_inst] = final_time;

    if(cb_max_times[last_inst] < final_time )
        cb_max_times[last_inst] = final_time;

    if( cb_operation ) {
     cb_times[last_inst] = (final_time +  times[last_inst]) / 2;
     cb_operation = false;
    } else {
     times[last_inst] = (final_time +  times[last_inst]) / 2; 
    }
}

void CPU::power2() {
  r.halt = false;
  r.stop = false;
  r.ei = false;
  r.ime = false;

  // normal instructions
  instructions[0x00] = &CPU::op_nop;
  instructions[0x01] = &CPU::op_ld_rr_nn<BC>;
  instructions[0x02] = &CPU::op_ld_rr_a<BC>;
  instructions[0x03] = &CPU::op_inc_rr<BC>;
  instructions[0x04] = &CPU::op_inc_r<B>;
  instructions[0x05] = &CPU::op_dec_r<B>;
  instructions[0x06] = &CPU::op_ld_r_n<B>;
  instructions[0x07] = &CPU::op_rlca;
  instructions[0x08] = &CPU::op_ld_nn_sp;
  instructions[0x09] = &CPU::op_add_hl_rr<BC>;
  instructions[0x0a] = &CPU::op_ld_a_rr<BC>;
  instructions[0x0b] = &CPU::op_dec_rr<BC>;
  instructions[0x0c] = &CPU::op_inc_r<C>;
  instructions[0x0d] = &CPU::op_dec_r<C>;
  instructions[0x0e] = &CPU::op_ld_r_n<C>;
  instructions[0x0f] = &CPU::op_rrca;
  instructions[0x10] = &CPU::op_stop;
  instructions[0x11] = &CPU::op_ld_rr_nn<DE>;
  instructions[0x12] = &CPU::op_ld_rr_a<DE>;
  instructions[0x13] = &CPU::op_inc_rr<DE>;
  instructions[0x14] = &CPU::op_inc_r<D>;
  instructions[0x15] = &CPU::op_dec_r<D>;
  instructions[0x16] = &CPU::op_ld_r_n<D>;
  instructions[0x17] = &CPU::op_rla;
  instructions[0x18] = &CPU::op_jr_n;
  instructions[0x19] = &CPU::op_add_hl_rr<DE>;
  instructions[0x1a] = &CPU::op_ld_a_rr<DE>;
  instructions[0x1b] = &CPU::op_dec_rr<DE>;
  instructions[0x1c] = &CPU::op_inc_r<E>;
  instructions[0x1d] = &CPU::op_dec_r<E>;
  instructions[0x1e] = &CPU::op_ld_r_n<E>;
  instructions[0x1f] = &CPU::op_rra;
  instructions[0x20] = &CPU::op_jr_f_n<ZF, 0>;
  instructions[0x21] = &CPU::op_ld_rr_nn<HL>;
  instructions[0x22] = &CPU::op_ldi_hl_a;
  instructions[0x23] = &CPU::op_inc_rr<HL>;
  instructions[0x24] = &CPU::op_inc_r<H>;
  instructions[0x25] = &CPU::op_dec_r<H>;
  instructions[0x26] = &CPU::op_ld_r_n<H>;
  instructions[0x27] = &CPU::op_daa;
  instructions[0x28] = &CPU::op_jr_f_n<ZF, 1>;
  instructions[0x29] = &CPU::op_add_hl_rr<HL>;
  instructions[0x2a] = &CPU::op_ldi_a_hl;
  instructions[0x2b] = &CPU::op_dec_rr<HL>;
  instructions[0x2c] = &CPU::op_inc_r<L>;
  instructions[0x2d] = &CPU::op_dec_r<L>;
  instructions[0x2e] = &CPU::op_ld_r_n<L>;
  instructions[0x2f] = &CPU::op_cpl;
  instructions[0x30] = &CPU::op_jr_f_n<CF, 0>;
  instructions[0x31] = &CPU::op_ld_rr_nn<SP>;
  instructions[0x32] = &CPU::op_ldd_hl_a;
  instructions[0x33] = &CPU::op_inc_rr<SP>;
  instructions[0x34] = &CPU::op_inc_hl;
  instructions[0x35] = &CPU::op_dec_hl;
  instructions[0x36] = &CPU::op_ld_hl_n;
  instructions[0x37] = &CPU::op_scf;
  instructions[0x38] = &CPU::op_jr_f_n<CF, 1>;
  instructions[0x39] = &CPU::op_add_hl_rr<SP>;
  instructions[0x3a] = &CPU::op_ldd_a_hl;
  instructions[0x3b] = &CPU::op_dec_rr<SP>;
  instructions[0x3c] = &CPU::op_inc_r<A>;
  instructions[0x3d] = &CPU::op_dec_r<A>;
  instructions[0x3e] = &CPU::op_ld_r_n<A>;
  instructions[0x3f] = &CPU::op_ccf;
  instructions[0x40] = &CPU::op_ld_r_r<B, B>;
  instructions[0x41] = &CPU::op_ld_r_r<B, C>;
  instructions[0x42] = &CPU::op_ld_r_r<B, D>;
  instructions[0x43] = &CPU::op_ld_r_r<B, E>;
  instructions[0x44] = &CPU::op_ld_r_r<B, H>;
  instructions[0x45] = &CPU::op_ld_r_r<B, L>;
  instructions[0x46] = &CPU::op_ld_r_hl<B>;
  instructions[0x47] = &CPU::op_ld_r_r<B, A>;
  instructions[0x48] = &CPU::op_ld_r_r<C, B>;
  instructions[0x49] = &CPU::op_ld_r_r<C, C>;
  instructions[0x4a] = &CPU::op_ld_r_r<C, D>;
  instructions[0x4b] = &CPU::op_ld_r_r<C, E>;
  instructions[0x4c] = &CPU::op_ld_r_r<C, H>;
  instructions[0x4d] = &CPU::op_ld_r_r<C, L>;
  instructions[0x4e] = &CPU::op_ld_r_hl<C>;
  instructions[0x4f] = &CPU::op_ld_r_r<C, A>;
  instructions[0x50] = &CPU::op_ld_r_r<D, B>;
  instructions[0x51] = &CPU::op_ld_r_r<D, C>;
  instructions[0x52] = &CPU::op_ld_r_r<D, D>;
  instructions[0x53] = &CPU::op_ld_r_r<D, E>;
  instructions[0x54] = &CPU::op_ld_r_r<D, H>;
  instructions[0x55] = &CPU::op_ld_r_r<D, L>;
  instructions[0x56] = &CPU::op_ld_r_hl<D>;
  instructions[0x57] = &CPU::op_ld_r_r<D, A>;
  instructions[0x58] = &CPU::op_ld_r_r<E, B>;
  instructions[0x59] = &CPU::op_ld_r_r<E, C>;
  instructions[0x5a] = &CPU::op_ld_r_r<E, D>;
  instructions[0x5b] = &CPU::op_ld_r_r<E, E>;
  instructions[0x5c] = &CPU::op_ld_r_r<E, H>;
  instructions[0x5d] = &CPU::op_ld_r_r<E, L>;
  instructions[0x5e] = &CPU::op_ld_r_hl<E>;
  instructions[0x5f] = &CPU::op_ld_r_r<E, A>;
  instructions[0x60] = &CPU::op_ld_r_r<H, B>;
  instructions[0x61] = &CPU::op_ld_r_r<H, C>;
  instructions[0x62] = &CPU::op_ld_r_r<H, D>;
  instructions[0x63] = &CPU::op_ld_r_r<H, E>;
  instructions[0x64] = &CPU::op_ld_r_r<H, H>;
  instructions[0x65] = &CPU::op_ld_r_r<H, L>;
  instructions[0x66] = &CPU::op_ld_r_hl<H>;
  instructions[0x67] = &CPU::op_ld_r_r<H, A>;
  instructions[0x68] = &CPU::op_ld_r_r<L, B>;
  instructions[0x69] = &CPU::op_ld_r_r<L, C>;
  instructions[0x6a] = &CPU::op_ld_r_r<L, D>;
  instructions[0x6b] = &CPU::op_ld_r_r<L, E>;
  instructions[0x6c] = &CPU::op_ld_r_r<L, H>;
  instructions[0x6d] = &CPU::op_ld_r_r<L, L>;
  instructions[0x6e] = &CPU::op_ld_r_hl<L>;
  instructions[0x6f] = &CPU::op_ld_r_r<L, A>;
  instructions[0x70] = &CPU::op_ld_hl_r<B>;
  instructions[0x71] = &CPU::op_ld_hl_r<C>;
  instructions[0x72] = &CPU::op_ld_hl_r<D>;
  instructions[0x73] = &CPU::op_ld_hl_r<E>;
  instructions[0x74] = &CPU::op_ld_hl_r<H>;
  instructions[0x75] = &CPU::op_ld_hl_r<L>;
  instructions[0x76] = &CPU::op_halt;
  instructions[0x77] = &CPU::op_ld_hl_r<A>;
  instructions[0x78] = &CPU::op_ld_r_r<A, B>;
  instructions[0x79] = &CPU::op_ld_r_r<A, C>;
  instructions[0x7a] = &CPU::op_ld_r_r<A, D>;
  instructions[0x7b] = &CPU::op_ld_r_r<A, E>;
  instructions[0x7c] = &CPU::op_ld_r_r<A, H>;
  instructions[0x7d] = &CPU::op_ld_r_r<A, L>;
  instructions[0x7e] = &CPU::op_ld_r_hl<A>;
  instructions[0x7f] = &CPU::op_ld_r_r<A, A>;
  instructions[0x80] = &CPU::op_add_a_r<B>;
  instructions[0x81] = &CPU::op_add_a_r<C>;
  instructions[0x82] = &CPU::op_add_a_r<D>;
  instructions[0x83] = &CPU::op_add_a_r<E>;
  instructions[0x84] = &CPU::op_add_a_r<H>;
  instructions[0x85] = &CPU::op_add_a_r<L>;
  instructions[0x86] = &CPU::op_add_a_hl;
  instructions[0x87] = &CPU::op_add_a_r<A>;
  instructions[0x88] = &CPU::op_adc_a_r<B>;
  instructions[0x89] = &CPU::op_adc_a_r<C>;
  instructions[0x8a] = &CPU::op_adc_a_r<D>;
  instructions[0x8b] = &CPU::op_adc_a_r<E>;
  instructions[0x8c] = &CPU::op_adc_a_r<H>;
  instructions[0x8d] = &CPU::op_adc_a_r<L>;
  instructions[0x8e] = &CPU::op_adc_a_hl;
  instructions[0x8f] = &CPU::op_adc_a_r<A>;
  instructions[0x90] = &CPU::op_sub_a_r<B>;
  instructions[0x91] = &CPU::op_sub_a_r<C>;
  instructions[0x92] = &CPU::op_sub_a_r<D>;
  instructions[0x93] = &CPU::op_sub_a_r<E>;
  instructions[0x94] = &CPU::op_sub_a_r<H>;
  instructions[0x95] = &CPU::op_sub_a_r<L>;
  instructions[0x96] = &CPU::op_sub_a_hl;
  instructions[0x97] = &CPU::op_sub_a_r<A>;
  instructions[0x98] = &CPU::op_sbc_a_r<B>;
  instructions[0x99] = &CPU::op_sbc_a_r<C>;
  instructions[0x9a] = &CPU::op_sbc_a_r<D>;
  instructions[0x9b] = &CPU::op_sbc_a_r<E>;
  instructions[0x9c] = &CPU::op_sbc_a_r<H>;
  instructions[0x9d] = &CPU::op_sbc_a_r<L>;
  instructions[0x9e] = &CPU::op_sbc_a_hl;
  instructions[0x9f] = &CPU::op_sbc_a_r<A>;
  instructions[0xa0] = &CPU::op_and_a_r<B>;
  instructions[0xa1] = &CPU::op_and_a_r<C>;
  instructions[0xa2] = &CPU::op_and_a_r<D>;
  instructions[0xa3] = &CPU::op_and_a_r<E>;
  instructions[0xa4] = &CPU::op_and_a_r<H>;
  instructions[0xa5] = &CPU::op_and_a_r<L>;
  instructions[0xa6] = &CPU::op_and_a_hl;
  instructions[0xa7] = &CPU::op_and_a_r<A>;
  instructions[0xa8] = &CPU::op_xor_a_r<B>;
  instructions[0xa9] = &CPU::op_xor_a_r<C>;
  instructions[0xaa] = &CPU::op_xor_a_r<D>;
  instructions[0xab] = &CPU::op_xor_a_r<E>;
  instructions[0xac] = &CPU::op_xor_a_r<H>;
  instructions[0xad] = &CPU::op_xor_a_r<L>;
  instructions[0xae] = &CPU::op_xor_a_hl;
  instructions[0xaf] = &CPU::op_xor_a_r<A>;
  instructions[0xb0] = &CPU::op_or_a_r<B>;
  instructions[0xb1] = &CPU::op_or_a_r<C>;
  instructions[0xb2] = &CPU::op_or_a_r<D>;
  instructions[0xb3] = &CPU::op_or_a_r<E>;
  instructions[0xb4] = &CPU::op_or_a_r<H>;
  instructions[0xb5] = &CPU::op_or_a_r<L>;
  instructions[0xb6] = &CPU::op_or_a_hl;
  instructions[0xb7] = &CPU::op_or_a_r<A>;
  instructions[0xb8] = &CPU::op_cp_a_r<B>;
  instructions[0xb9] = &CPU::op_cp_a_r<C>;
  instructions[0xba] = &CPU::op_cp_a_r<D>;
  instructions[0xbb] = &CPU::op_cp_a_r<E>;
  instructions[0xbc] = &CPU::op_cp_a_r<H>;
  instructions[0xbd] = &CPU::op_cp_a_r<L>;
  instructions[0xbe] = &CPU::op_cp_a_hl;
  instructions[0xbf] = &CPU::op_cp_a_r<A>;
  instructions[0xc0] = &CPU::op_ret_f<ZF, 0>;
  instructions[0xc1] = &CPU::op_pop_rr<BC>;
  instructions[0xc2] = &CPU::op_jp_f_nn<ZF, 0>;
  instructions[0xc3] = &CPU::op_jp_nn;
  instructions[0xc4] = &CPU::op_call_f_nn<ZF, 0>;
  instructions[0xc5] = &CPU::op_push_rr<BC>;
  instructions[0xc6] = &CPU::op_add_a_n;
  instructions[0xc7] = &CPU::op_rst_n<0x00>;
  instructions[0xc8] = &CPU::op_ret_f<ZF, 1>;
  instructions[0xc9] = &CPU::op_ret;
  instructions[0xca] = &CPU::op_jp_f_nn<ZF, 1>;
  instructions[0xcb] = &CPU::op_cb;
  instructions[0xcc] = &CPU::op_call_f_nn<ZF, 1>;
  instructions[0xcd] = &CPU::op_call_nn;
  instructions[0xce] = &CPU::op_adc_a_n;
  instructions[0xcf] = &CPU::op_rst_n<0x08>;
  instructions[0xd0] = &CPU::op_ret_f<CF, 0>;
  instructions[0xd1] = &CPU::op_pop_rr<DE>;
  instructions[0xd2] = &CPU::op_jp_f_nn<CF, 0>;
  instructions[0xd3] = &CPU::op_xx;
  instructions[0xd4] = &CPU::op_call_f_nn<CF, 0>;
  instructions[0xd5] = &CPU::op_push_rr<DE>;
  instructions[0xd6] = &CPU::op_sub_a_n;
  instructions[0xd7] = &CPU::op_rst_n<0x10>;
  instructions[0xd8] = &CPU::op_ret_f<CF, 1>;
  instructions[0xd9] = &CPU::op_reti;
  instructions[0xda] = &CPU::op_jp_f_nn<CF, 1>;
  instructions[0xdb] = &CPU::op_xx;
  instructions[0xdc] = &CPU::op_call_f_nn<CF, 1>;
  instructions[0xdd] = &CPU::op_xx;
  instructions[0xde] = &CPU::op_sbc_a_n;
  instructions[0xdf] = &CPU::op_rst_n<0x18>;
  instructions[0xe0] = &CPU::op_ld_ffn_a;
  instructions[0xe1] = &CPU::op_pop_rr<HL>;
  instructions[0xe2] = &CPU::op_ld_ffc_a;
  instructions[0xe3] = &CPU::op_xx;
  instructions[0xe4] = &CPU::op_xx;
  instructions[0xe5] = &CPU::op_push_rr<HL>;
  instructions[0xe6] = &CPU::op_and_a_n;
  instructions[0xe7] = &CPU::op_rst_n<0x20>;
  instructions[0xe8] = &CPU::op_add_sp_n;
  instructions[0xe9] = &CPU::op_jp_hl;
  instructions[0xea] = &CPU::op_ld_nn_a;
  instructions[0xeb] = &CPU::op_xx;
  instructions[0xec] = &CPU::op_xx;
  instructions[0xed] = &CPU::op_xx;
  instructions[0xee] = &CPU::op_xor_a_n;
  instructions[0xef] = &CPU::op_rst_n<0x28>;
  instructions[0xf0] = &CPU::op_ld_a_ffn;
  instructions[0xf1] = &CPU::op_pop_rr<AF>;
  instructions[0xf2] = &CPU::op_ld_a_ffc;
  instructions[0xf3] = &CPU::op_di;
  instructions[0xf4] = &CPU::op_xx;
  instructions[0xf5] = &CPU::op_push_rr<AF>;
  instructions[0xf6] = &CPU::op_or_a_n;
  instructions[0xf7] = &CPU::op_rst_n<0x30>;
  instructions[0xf8] = &CPU::op_ld_hl_sp_n;
  instructions[0xf9] = &CPU::op_ld_sp_hl;
  instructions[0xfa] = &CPU::op_ld_a_nn;
  instructions[0xfb] = &CPU::op_ei;
  instructions[0xfc] = &CPU::op_xx;
  instructions[0xfd] = &CPU::op_xx;
  instructions[0xfe] = &CPU::op_cp_a_n;
  instructions[0xff] = &CPU::op_rst_n<0x38>;

  // cb instructions
  cb_instructions[0x00] = &CPU::op_rlc_r<B>;
  cb_instructions[0x01] = &CPU::op_rlc_r<C>;
  cb_instructions[0x02] = &CPU::op_rlc_r<D>;
  cb_instructions[0x03] = &CPU::op_rlc_r<E>;
  cb_instructions[0x04] = &CPU::op_rlc_r<H>;
  cb_instructions[0x05] = &CPU::op_rlc_r<L>;
  cb_instructions[0x06] = &CPU::op_rlc_hl;
  cb_instructions[0x07] = &CPU::op_rlc_r<A>;
  cb_instructions[0x08] = &CPU::op_rrc_r<B>;
  cb_instructions[0x09] = &CPU::op_rrc_r<C>;
  cb_instructions[0x0a] = &CPU::op_rrc_r<D>;
  cb_instructions[0x0b] = &CPU::op_rrc_r<E>;
  cb_instructions[0x0c] = &CPU::op_rrc_r<H>;
  cb_instructions[0x0d] = &CPU::op_rrc_r<L>;
  cb_instructions[0x0e] = &CPU::op_rrc_hl;
  cb_instructions[0x0f] = &CPU::op_rrc_r<A>;
  cb_instructions[0x10] = &CPU::op_rl_r<B>;
  cb_instructions[0x11] = &CPU::op_rl_r<C>;
  cb_instructions[0x12] = &CPU::op_rl_r<D>;
  cb_instructions[0x13] = &CPU::op_rl_r<E>;
  cb_instructions[0x14] = &CPU::op_rl_r<H>;
  cb_instructions[0x15] = &CPU::op_rl_r<L>;
  cb_instructions[0x16] = &CPU::op_rl_hl;
  cb_instructions[0x17] = &CPU::op_rl_r<A>;
  cb_instructions[0x18] = &CPU::op_rr_r<B>;
  cb_instructions[0x19] = &CPU::op_rr_r<C>;
  cb_instructions[0x1a] = &CPU::op_rr_r<D>;
  cb_instructions[0x1b] = &CPU::op_rr_r<E>;
  cb_instructions[0x1c] = &CPU::op_rr_r<H>;
  cb_instructions[0x1d] = &CPU::op_rr_r<L>;
  cb_instructions[0x1e] = &CPU::op_rr_hl;
  cb_instructions[0x1f] = &CPU::op_rr_r<A>;
  cb_instructions[0x20] = &CPU::op_sla_r<B>;
  cb_instructions[0x21] = &CPU::op_sla_r<C>;
  cb_instructions[0x22] = &CPU::op_sla_r<D>;
  cb_instructions[0x23] = &CPU::op_sla_r<E>;
  cb_instructions[0x24] = &CPU::op_sla_r<H>;
  cb_instructions[0x25] = &CPU::op_sla_r<L>;
  cb_instructions[0x26] = &CPU::op_sla_hl;
  cb_instructions[0x27] = &CPU::op_sla_r<A>;
  cb_instructions[0x28] = &CPU::op_sra_r<B>;
  cb_instructions[0x29] = &CPU::op_sra_r<C>;
  cb_instructions[0x2a] = &CPU::op_sra_r<D>;
  cb_instructions[0x2b] = &CPU::op_sra_r<E>;
  cb_instructions[0x2c] = &CPU::op_sra_r<H>;
  cb_instructions[0x2d] = &CPU::op_sra_r<L>;
  cb_instructions[0x2e] = &CPU::op_sra_hl;
  cb_instructions[0x2f] = &CPU::op_sra_r<A>;
  cb_instructions[0x30] = &CPU::op_swap_r<B>;
  cb_instructions[0x31] = &CPU::op_swap_r<C>;
  cb_instructions[0x32] = &CPU::op_swap_r<D>;
  cb_instructions[0x33] = &CPU::op_swap_r<E>;
  cb_instructions[0x34] = &CPU::op_swap_r<H>;
  cb_instructions[0x35] = &CPU::op_swap_r<L>;
  cb_instructions[0x36] = &CPU::op_swap_hl;
  cb_instructions[0x37] = &CPU::op_swap_r<A>;
  cb_instructions[0x38] = &CPU::op_srl_r<B>;
  cb_instructions[0x39] = &CPU::op_srl_r<C>;
  cb_instructions[0x3a] = &CPU::op_srl_r<D>;
  cb_instructions[0x3b] = &CPU::op_srl_r<E>;
  cb_instructions[0x3c] = &CPU::op_srl_r<H>;
  cb_instructions[0x3d] = &CPU::op_srl_r<L>;
  cb_instructions[0x3e] = &CPU::op_srl_hl;
  cb_instructions[0x3f] = &CPU::op_srl_r<A>;
  cb_instructions[0x40] = &CPU::op_bit_n_r<0, B>;
  cb_instructions[0x41] = &CPU::op_bit_n_r<0, C>;
  cb_instructions[0x42] = &CPU::op_bit_n_r<0, D>;
  cb_instructions[0x43] = &CPU::op_bit_n_r<0, E>;
  cb_instructions[0x44] = &CPU::op_bit_n_r<0, H>;
  cb_instructions[0x45] = &CPU::op_bit_n_r<0, L>;
  cb_instructions[0x46] = &CPU::op_bit_n_hl<0>;
  cb_instructions[0x47] = &CPU::op_bit_n_r<0, A>;
  cb_instructions[0x48] = &CPU::op_bit_n_r<1, B>;
  cb_instructions[0x49] = &CPU::op_bit_n_r<1, C>;
  cb_instructions[0x4a] = &CPU::op_bit_n_r<1, D>;
  cb_instructions[0x4b] = &CPU::op_bit_n_r<1, E>;
  cb_instructions[0x4c] = &CPU::op_bit_n_r<1, H>;
  cb_instructions[0x4d] = &CPU::op_bit_n_r<1, L>;
  cb_instructions[0x4e] = &CPU::op_bit_n_hl<1>;
  cb_instructions[0x4f] = &CPU::op_bit_n_r<1, A>;
  cb_instructions[0x50] = &CPU::op_bit_n_r<2, B>;
  cb_instructions[0x51] = &CPU::op_bit_n_r<2, C>;
  cb_instructions[0x52] = &CPU::op_bit_n_r<2, D>;
  cb_instructions[0x53] = &CPU::op_bit_n_r<2, E>;
  cb_instructions[0x54] = &CPU::op_bit_n_r<2, H>;
  cb_instructions[0x55] = &CPU::op_bit_n_r<2, L>;
  cb_instructions[0x56] = &CPU::op_bit_n_hl<2>;
  cb_instructions[0x57] = &CPU::op_bit_n_r<2, A>;
  cb_instructions[0x58] = &CPU::op_bit_n_r<3, B>;
  cb_instructions[0x59] = &CPU::op_bit_n_r<3, C>;
  cb_instructions[0x5a] = &CPU::op_bit_n_r<3, D>;
  cb_instructions[0x5b] = &CPU::op_bit_n_r<3, E>;
  cb_instructions[0x5c] = &CPU::op_bit_n_r<3, H>;
  cb_instructions[0x5d] = &CPU::op_bit_n_r<3, L>;
  cb_instructions[0x5e] = &CPU::op_bit_n_hl<3>;
  cb_instructions[0x5f] = &CPU::op_bit_n_r<3, A>;
  cb_instructions[0x60] = &CPU::op_bit_n_r<4, B>;
  cb_instructions[0x61] = &CPU::op_bit_n_r<4, C>;
  cb_instructions[0x62] = &CPU::op_bit_n_r<4, D>;
  cb_instructions[0x63] = &CPU::op_bit_n_r<4, E>;
  cb_instructions[0x64] = &CPU::op_bit_n_r<4, H>;
  cb_instructions[0x65] = &CPU::op_bit_n_r<4, L>;
  cb_instructions[0x66] = &CPU::op_bit_n_hl<4>;
  cb_instructions[0x67] = &CPU::op_bit_n_r<4, A>;
  cb_instructions[0x68] = &CPU::op_bit_n_r<5, B>;
  cb_instructions[0x69] = &CPU::op_bit_n_r<5, C>;
  cb_instructions[0x6a] = &CPU::op_bit_n_r<5, D>;
  cb_instructions[0x6b] = &CPU::op_bit_n_r<5, E>;
  cb_instructions[0x6c] = &CPU::op_bit_n_r<5, H>;
  cb_instructions[0x6d] = &CPU::op_bit_n_r<5, L>;
  cb_instructions[0x6e] = &CPU::op_bit_n_hl<5>;
  cb_instructions[0x6f] = &CPU::op_bit_n_r<5, A>;
  cb_instructions[0x70] = &CPU::op_bit_n_r<6, B>;
  cb_instructions[0x71] = &CPU::op_bit_n_r<6, C>;
  cb_instructions[0x72] = &CPU::op_bit_n_r<6, D>;
  cb_instructions[0x73] = &CPU::op_bit_n_r<6, E>;
  cb_instructions[0x74] = &CPU::op_bit_n_r<6, H>;
  cb_instructions[0x75] = &CPU::op_bit_n_r<6, L>;
  cb_instructions[0x76] = &CPU::op_bit_n_hl<6>;
  cb_instructions[0x77] = &CPU::op_bit_n_r<6, A>;
  cb_instructions[0x78] = &CPU::op_bit_n_r<7, B>;
  cb_instructions[0x79] = &CPU::op_bit_n_r<7, C>;
  cb_instructions[0x7a] = &CPU::op_bit_n_r<7, D>;
  cb_instructions[0x7b] = &CPU::op_bit_n_r<7, E>;
  cb_instructions[0x7c] = &CPU::op_bit_n_r<7, H>;
  cb_instructions[0x7d] = &CPU::op_bit_n_r<7, L>;
  cb_instructions[0x7e] = &CPU::op_bit_n_hl<7>;
  cb_instructions[0x7f] = &CPU::op_bit_n_r<7, A>;
  cb_instructions[0x80] = &CPU::op_res_n_r<0, B>;
  cb_instructions[0x81] = &CPU::op_res_n_r<0, C>;
  cb_instructions[0x82] = &CPU::op_res_n_r<0, D>;
  cb_instructions[0x83] = &CPU::op_res_n_r<0, E>;
  cb_instructions[0x84] = &CPU::op_res_n_r<0, H>;
  cb_instructions[0x85] = &CPU::op_res_n_r<0, L>;
  cb_instructions[0x86] = &CPU::op_res_n_hl<0>;
  cb_instructions[0x87] = &CPU::op_res_n_r<0, A>;
  cb_instructions[0x88] = &CPU::op_res_n_r<1, B>;
  cb_instructions[0x89] = &CPU::op_res_n_r<1, C>;
  cb_instructions[0x8a] = &CPU::op_res_n_r<1, D>;
  cb_instructions[0x8b] = &CPU::op_res_n_r<1, E>;
  cb_instructions[0x8c] = &CPU::op_res_n_r<1, H>;
  cb_instructions[0x8d] = &CPU::op_res_n_r<1, L>;
  cb_instructions[0x8e] = &CPU::op_res_n_hl<1>;
  cb_instructions[0x8f] = &CPU::op_res_n_r<1, A>;
  cb_instructions[0x90] = &CPU::op_res_n_r<2, B>;
  cb_instructions[0x91] = &CPU::op_res_n_r<2, C>;
  cb_instructions[0x92] = &CPU::op_res_n_r<2, D>;
  cb_instructions[0x93] = &CPU::op_res_n_r<2, E>;
  cb_instructions[0x94] = &CPU::op_res_n_r<2, H>;
  cb_instructions[0x95] = &CPU::op_res_n_r<2, L>;
  cb_instructions[0x96] = &CPU::op_res_n_hl<2>;
  cb_instructions[0x97] = &CPU::op_res_n_r<2, A>;
  cb_instructions[0x98] = &CPU::op_res_n_r<3, B>;
  cb_instructions[0x99] = &CPU::op_res_n_r<3, C>;
  cb_instructions[0x9a] = &CPU::op_res_n_r<3, D>;
  cb_instructions[0x9b] = &CPU::op_res_n_r<3, E>;
  cb_instructions[0x9c] = &CPU::op_res_n_r<3, H>;
  cb_instructions[0x9d] = &CPU::op_res_n_r<3, L>;
  cb_instructions[0x9e] = &CPU::op_res_n_hl<3>;
  cb_instructions[0x9f] = &CPU::op_res_n_r<3, A>;
  cb_instructions[0xa0] = &CPU::op_res_n_r<4, B>;
  cb_instructions[0xa1] = &CPU::op_res_n_r<4, C>;
  cb_instructions[0xa2] = &CPU::op_res_n_r<4, D>;
  cb_instructions[0xa3] = &CPU::op_res_n_r<4, E>;
  cb_instructions[0xa4] = &CPU::op_res_n_r<4, H>;
  cb_instructions[0xa5] = &CPU::op_res_n_r<4, L>;
  cb_instructions[0xa6] = &CPU::op_res_n_hl<4>;
  cb_instructions[0xa7] = &CPU::op_res_n_r<4, A>;
  cb_instructions[0xa8] = &CPU::op_res_n_r<5, B>;
  cb_instructions[0xa9] = &CPU::op_res_n_r<5, C>;
  cb_instructions[0xaa] = &CPU::op_res_n_r<5, D>;
  cb_instructions[0xab] = &CPU::op_res_n_r<5, E>;
  cb_instructions[0xac] = &CPU::op_res_n_r<5, H>;
  cb_instructions[0xad] = &CPU::op_res_n_r<5, L>;
  cb_instructions[0xae] = &CPU::op_res_n_hl<5>;
  cb_instructions[0xaf] = &CPU::op_res_n_r<5, A>;
  cb_instructions[0xb0] = &CPU::op_res_n_r<6, B>;
  cb_instructions[0xb1] = &CPU::op_res_n_r<6, C>;
  cb_instructions[0xb2] = &CPU::op_res_n_r<6, D>;
  cb_instructions[0xb3] = &CPU::op_res_n_r<6, E>;
  cb_instructions[0xb4] = &CPU::op_res_n_r<6, H>;
  cb_instructions[0xb5] = &CPU::op_res_n_r<6, L>;
  cb_instructions[0xb6] = &CPU::op_res_n_hl<6>;
  cb_instructions[0xb7] = &CPU::op_res_n_r<6, A>;
  cb_instructions[0xb8] = &CPU::op_res_n_r<7, B>;
  cb_instructions[0xb9] = &CPU::op_res_n_r<7, C>;
  cb_instructions[0xba] = &CPU::op_res_n_r<7, D>;
  cb_instructions[0xbb] = &CPU::op_res_n_r<7, E>;
  cb_instructions[0xbc] = &CPU::op_res_n_r<7, H>;
  cb_instructions[0xbd] = &CPU::op_res_n_r<7, L>;
  cb_instructions[0xbe] = &CPU::op_res_n_hl<7>;
  cb_instructions[0xbf] = &CPU::op_res_n_r<7, A>;
  cb_instructions[0xc0] = &CPU::op_set_n_r<0, B>;
  cb_instructions[0xc1] = &CPU::op_set_n_r<0, C>;
  cb_instructions[0xc2] = &CPU::op_set_n_r<0, D>;
  cb_instructions[0xc3] = &CPU::op_set_n_r<0, E>;
  cb_instructions[0xc4] = &CPU::op_set_n_r<0, H>;
  cb_instructions[0xc5] = &CPU::op_set_n_r<0, L>;
  cb_instructions[0xc6] = &CPU::op_set_n_hl<0>;
  cb_instructions[0xc7] = &CPU::op_set_n_r<0, A>;
  cb_instructions[0xc8] = &CPU::op_set_n_r<1, B>;
  cb_instructions[0xc9] = &CPU::op_set_n_r<1, C>;
  cb_instructions[0xca] = &CPU::op_set_n_r<1, D>;
  cb_instructions[0xcb] = &CPU::op_set_n_r<1, E>;
  cb_instructions[0xcc] = &CPU::op_set_n_r<1, H>;
  cb_instructions[0xcd] = &CPU::op_set_n_r<1, L>;
  cb_instructions[0xce] = &CPU::op_set_n_hl<1>;
  cb_instructions[0xcf] = &CPU::op_set_n_r<1, A>;
  cb_instructions[0xd0] = &CPU::op_set_n_r<2, B>;
  cb_instructions[0xd1] = &CPU::op_set_n_r<2, C>;
  cb_instructions[0xd2] = &CPU::op_set_n_r<2, D>;
  cb_instructions[0xd3] = &CPU::op_set_n_r<2, E>;
  cb_instructions[0xd4] = &CPU::op_set_n_r<2, H>;
  cb_instructions[0xd5] = &CPU::op_set_n_r<2, L>;
  cb_instructions[0xd6] = &CPU::op_set_n_hl<2>;
  cb_instructions[0xd7] = &CPU::op_set_n_r<2, A>;
  cb_instructions[0xd8] = &CPU::op_set_n_r<3, B>;
  cb_instructions[0xd9] = &CPU::op_set_n_r<3, C>;
  cb_instructions[0xda] = &CPU::op_set_n_r<3, D>;
  cb_instructions[0xdb] = &CPU::op_set_n_r<3, E>;
  cb_instructions[0xdc] = &CPU::op_set_n_r<3, H>;
  cb_instructions[0xdd] = &CPU::op_set_n_r<3, L>;
  cb_instructions[0xde] = &CPU::op_set_n_hl<3>;
  cb_instructions[0xdf] = &CPU::op_set_n_r<3, A>;
  cb_instructions[0xe0] = &CPU::op_set_n_r<4, B>;
  cb_instructions[0xe1] = &CPU::op_set_n_r<4, C>;
  cb_instructions[0xe2] = &CPU::op_set_n_r<4, D>;
  cb_instructions[0xe3] = &CPU::op_set_n_r<4, E>;
  cb_instructions[0xe4] = &CPU::op_set_n_r<4, H>;
  cb_instructions[0xe5] = &CPU::op_set_n_r<4, L>;
  cb_instructions[0xe6] = &CPU::op_set_n_hl<4>;
  cb_instructions[0xe7] = &CPU::op_set_n_r<4, A>;
  cb_instructions[0xe8] = &CPU::op_set_n_r<5, B>;
  cb_instructions[0xe9] = &CPU::op_set_n_r<5, C>;
  cb_instructions[0xea] = &CPU::op_set_n_r<5, D>;
  cb_instructions[0xeb] = &CPU::op_set_n_r<5, E>;
  cb_instructions[0xec] = &CPU::op_set_n_r<5, H>;
  cb_instructions[0xed] = &CPU::op_set_n_r<5, L>;
  cb_instructions[0xee] = &CPU::op_set_n_hl<5>;
  cb_instructions[0xef] = &CPU::op_set_n_r<5, A>;
  cb_instructions[0xf0] = &CPU::op_set_n_r<6, B>;
  cb_instructions[0xf1] = &CPU::op_set_n_r<6, C>;
  cb_instructions[0xf2] = &CPU::op_set_n_r<6, D>;
  cb_instructions[0xf3] = &CPU::op_set_n_r<6, E>;
  cb_instructions[0xf4] = &CPU::op_set_n_r<6, H>;
  cb_instructions[0xf5] = &CPU::op_set_n_r<6, L>;
  cb_instructions[0xf6] = &CPU::op_set_n_hl<6>;
  cb_instructions[0xf7] = &CPU::op_set_n_r<6, A>;
  cb_instructions[0xf8] = &CPU::op_set_n_r<7, B>;
  cb_instructions[0xf9] = &CPU::op_set_n_r<7, C>;
  cb_instructions[0xfa] = &CPU::op_set_n_r<7, D>;
  cb_instructions[0xfb] = &CPU::op_set_n_r<7, E>;
  cb_instructions[0xfc] = &CPU::op_set_n_r<7, H>;
  cb_instructions[0xfd] = &CPU::op_set_n_r<7, L>;
  cb_instructions[0xfe] = &CPU::op_set_n_hl<7>;
  cb_instructions[0xff] = &CPU::op_set_n_r<7, A>;
}

void CPU::exec_cb() {
  cb_operation = true;
  uint8 opcode = op_read(r[PC]++);
  last_inst = opcode;
  cb_inst_counter[opcode]++;
  (this->*cb_instructions[opcode])();
}

void CPU::exec() {
  uint8 opcode = op_read(r[PC]++);
  last_inst = opcode;
  inst_counter[opcode]++;
  instruction_count++;
  (this->*instructions[opcode])();

}

void CPU::Main() {
  cpu.main();
}

void CPU::main() {
	std::clock_t time;

  while(true) {
    if(scheduler.sync == Scheduler::SynchronizeMode::CPU) {
      scheduler.sync = Scheduler::SynchronizeMode::All;
      scheduler.exit(Scheduler::ExitReason::SynchronizeEvent);
    }

    interrupt_test();

    time = std::clock();
    exec();
    double final_time = double(std::clock() - time) / (double)CLOCKS_PER_SEC;
    get_times(final_time);
  }
  
}

void CPU::interrupt_raise(CPU::Interrupt id) {
  if(id == Interrupt::Vblank) {
    status.interrupt_request_vblank = 1;
    if(status.interrupt_enable_vblank) r.halt = false;
  }

  if(id == Interrupt::Stat) {
    status.interrupt_request_stat = 1;
    if(status.interrupt_enable_stat) r.halt = false;
  }

  if(id == Interrupt::Timer) {
    status.interrupt_request_timer = 1;
    if(status.interrupt_enable_timer) r.halt = false;
  }

  if(id == Interrupt::Serial) {
    status.interrupt_request_serial = 1;
    if(status.interrupt_enable_serial) r.halt = false;
  }

  if(id == Interrupt::Joypad) {
    status.interrupt_request_joypad = 1;
    if(status.interrupt_enable_joypad) r.halt = r.stop = false;
  }
}

void CPU::interrupt_test() {
  if(r.ime) {
    if(status.interrupt_request_vblank && status.interrupt_enable_vblank) {
      status.interrupt_request_vblank = 0;
      return interrupt_exec(0x0040);
    }

    if(status.interrupt_request_stat && status.interrupt_enable_stat) {
      status.interrupt_request_stat = 0;
      return interrupt_exec(0x0048);
    }

    if(status.interrupt_request_timer && status.interrupt_enable_timer) {
      status.interrupt_request_timer = 0;
      return interrupt_exec(0x0050);
    }

    if(status.interrupt_request_serial && status.interrupt_enable_serial) {
      status.interrupt_request_serial = 0;
      return interrupt_exec(0x0058);
    }

    if(status.interrupt_request_joypad && status.interrupt_enable_joypad) {
      status.interrupt_request_joypad = 0;
      return interrupt_exec(0x0060);
    }
  }
}

void CPU::interrupt_exec(uint16 pc) {
  r.ime = 0;
  op_write(--r[SP], r[PC] >> 8);
  op_write(--r[SP], r[PC] >> 0);
  r[PC] = pc;
  op_io();
  op_io();
  op_io();
}

bool CPU::stop() {
  if(status.speed_switch) {
    status.speed_switch = 0;
    status.speed_double ^= 1;
    if(status.speed_double == 0) frequency = 4 * 1024 * 1024;
    if(status.speed_double == 1) frequency = 8 * 1024 * 1024;
    return true;
  }
  return false;
}

void CPU::power() {
  create(Main, 4 * 1024 * 1024);
  power2();

  for(unsigned n = 0xc000; n <= 0xdfff; n++) bus.mmio[n] = this;  //WRAM
  for(unsigned n = 0xe000; n <= 0xfdff; n++) bus.mmio[n] = this;  //WRAM (mirror)
  for(unsigned n = 0xff80; n <= 0xfffe; n++) bus.mmio[n] = this;  //HRAM

  bus.mmio[0xff00] = this;  //JOYP
  bus.mmio[0xff01] = this;  //SB
  bus.mmio[0xff02] = this;  //SC
  bus.mmio[0xff04] = this;  //DIV
  bus.mmio[0xff05] = this;  //TIMA
  bus.mmio[0xff06] = this;  //TMA
  bus.mmio[0xff07] = this;  //TAC
  bus.mmio[0xff0f] = this;  //IF
  bus.mmio[0xff46] = this;  //DMA
  bus.mmio[0xffff] = this;  //IE

  if(system.cgb()) {
  bus.mmio[0xff4d] = this;  //KEY1
  bus.mmio[0xff51] = this;  //HDMA1
  bus.mmio[0xff52] = this;  //HDMA2
  bus.mmio[0xff53] = this;  //HDMA3
  bus.mmio[0xff54] = this;  //HDMA4
  bus.mmio[0xff55] = this;  //HDMA5
  bus.mmio[0xff56] = this;  //RP
  bus.mmio[0xff6c] = this;  //???
  bus.mmio[0xff70] = this;  //SVBK
  bus.mmio[0xff72] = this;  //???
  bus.mmio[0xff73] = this;  //???
  bus.mmio[0xff74] = this;  //???
  bus.mmio[0xff75] = this;  //???
  bus.mmio[0xff76] = this;  //???
  bus.mmio[0xff77] = this;  //???
  }

  for(auto& n : wram) n = 0x00;
  for(auto& n : hram) n = 0x00;

  r[PC] = 0x0000;
  r[SP] = 0x0000;
  r[AF] = 0x0000;
  r[BC] = 0x0000;
  r[DE] = 0x0000;
  r[HL] = 0x0000;

  status.clock = 0;

  status.p15 = 0;
  status.p14 = 0;
  status.joyp = 0;
  status.mlt_req = 0;

  status.serial_data = 0;
  status.serial_bits = 0;

  status.serial_transfer = 0;
  status.serial_clock = 0;

  status.div = 0;

  status.tima = 0;

  status.tma = 0;

  status.timer_enable = 0;
  status.timer_clock = 0;

  status.interrupt_request_joypad = 0;
  status.interrupt_request_serial = 0;
  status.interrupt_request_timer = 0;
  status.interrupt_request_stat = 0;
  status.interrupt_request_vblank = 0;

  status.speed_double = 0;
  status.speed_switch = 0;

  status.dma_source = 0;
  status.dma_target = 0;

  status.dma_mode = 0;
  status.dma_length = 0;

  status.ff6c = 0;
  status.ff72 = 0;
  status.ff73 = 0;
  status.ff74 = 0;
  status.ff75 = 0;

  status.wram_bank = 1;

  status.interrupt_enable_joypad = 0;
  status.interrupt_enable_serial = 0;
  status.interrupt_enable_timer = 0;
  status.interrupt_enable_stat = 0;
  status.interrupt_enable_vblank = 0;
}

}

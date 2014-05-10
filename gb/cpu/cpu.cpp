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

// macro definitions
#define test_for_interrupt() if(r.ime){interrupt_test();} 

inline void CPU::op_xx() {
}

inline void CPU::op_cb() {
  exec_cb();
}

//8-bit load commands

template<unsigned x, unsigned y> inline void CPU::op_ld_r_r() {
  r[x] = r[y];
}

template<unsigned x> inline void CPU::op_ld_r_n() {
  r[x] = op_read(r[PC]++);
}

template<unsigned x> inline void CPU::op_ld_r_hl() {
  r[x] = op_read(r[HL]);
}

template<unsigned x> inline void CPU::op_ld_hl_r() {
  op_write(r[HL], r[x]);
}

inline void CPU::op_ld_hl_n() {
  op_write(r[HL], op_read(r[PC]++));
}

template<unsigned x> inline void CPU::op_ld_a_rr() {
  r[A] = op_read(r[x]);
}

inline void CPU::op_ld_a_nn() {
  uint8 lo = op_read(r[PC]++);
  uint8 hi = op_read(r[PC]++);
  r[A] = op_read((hi << 8) | (lo << 0));
}

template<unsigned x> inline void CPU::op_ld_rr_a() {
  op_write(r[x], r[A]);
}

inline void CPU::op_ld_nn_a() {
  uint8 lo = op_read(r[PC]++);
  uint8 hi = op_read(r[PC]++);
  op_write((hi << 8) | (lo << 0), r[A]);
}

inline void CPU::op_ld_a_ffn() {
  r[A] = op_read(0xff00 + op_read(r[PC]++));
}

inline void CPU::op_ld_ffn_a() {
  op_write(0xff00 + op_read(r[PC]++), r[A]);
}

inline void CPU::op_ld_a_ffc() {
  r[A] = op_read(0xff00 + r[C]);
}

inline void CPU::op_ld_ffc_a() {
  op_write(0xff00 + r[C], r[A]);
}

inline void CPU::op_ldi_hl_a() {
  op_write(r[HL], r[A]);
  r[HL]++;
}

inline void CPU::op_ldi_a_hl() {
  r[A] = op_read(r[HL]);
  r[HL]++;
}

inline void CPU::op_ldd_hl_a() {
  op_write(r[HL], r[A]);
  r[HL]--;
}

inline void CPU::op_ldd_a_hl() {
  r[A] = op_read(r[HL]);
  r[HL]--;
}

//16-bit load commands

template<unsigned x> inline void CPU::op_ld_rr_nn() {
  r[x]  = op_read(r[PC]++) << 0;
  r[x] |= op_read(r[PC]++) << 8;
}

inline void CPU::op_ld_nn_sp() {
  uint16 addr = op_read(r[PC]++) << 0;
  addr |= op_read(r[PC]++) << 8;
  op_write(addr + 0, r[SP] >> 0);
  op_write(addr + 1, r[SP] >> 8);
}

inline void CPU::op_ld_sp_hl() {
  r[SP] = r[HL];
  op_io();
}

template<unsigned x> inline void CPU::op_push_rr() {
  op_write(--r[SP], r[x] >> 8);
  op_write(--r[SP], r[x] >> 0);
  op_io();
}

template<unsigned x> inline void CPU::op_pop_rr() {
  r[x]  = op_read(r[SP]++) << 0;
  r[x] |= op_read(r[SP]++) << 8;
}

//8-bit arithmetic commands

inline void CPU::opi_add_a(uint8 x) {
  uint16 rh = r[A] + x;
  uint16 rl = (r[A] & 0x0f) + (x & 0x0f);
  r[A] = rh;
  r.f.z = (uint8)rh == 0;
  r.f.n = 0;
  r.f.h = rl > 0x0f;
  r.f.c = rh > 0xff;
}

template<unsigned x> inline void CPU::op_add_a_r() { opi_add_a(r[x]); }
inline void CPU::op_add_a_n() { opi_add_a(op_read(r[PC]++)); }
inline void CPU::op_add_a_hl() { opi_add_a(op_read(r[HL])); }

inline void CPU::opi_adc_a(uint8 x) {
  uint16 rh = r[A] + x + r.f.c;
  uint16 rl = (r[A] & 0x0f) + (x & 0x0f) + r.f.c;
  r[A] = rh;
  r.f.z = (uint8)rh == 0;
  r.f.n = 0;
  r.f.h = rl > 0x0f;
  r.f.c = rh > 0xff;
}

template<unsigned x> inline void CPU::op_adc_a_r() { opi_adc_a(r[x]); }
inline void CPU::op_adc_a_n() { opi_adc_a(op_read(r[PC]++)); }
inline void CPU::op_adc_a_hl() { opi_adc_a(op_read(r[HL])); }

inline void CPU::opi_sub_a(uint8 x) {
  uint16 rh = r[A] - x;
  uint16 rl = (r[A] & 0x0f) - (x & 0x0f);
  r[A] = rh;
  r.f.z = (uint8)rh == 0;
  r.f.n = 1;
  r.f.h = rl > 0x0f;
  r.f.c = rh > 0xff;
}

template<unsigned x> inline void CPU::op_sub_a_r() { opi_sub_a(r[x]); }
inline void CPU::op_sub_a_n() { opi_sub_a(op_read(r[PC]++)); }
inline void CPU::op_sub_a_hl() { opi_sub_a(op_read(r[HL])); }

inline void CPU::opi_sbc_a(uint8 x) {
  uint16 rh = r[A] - x - r.f.c;
  uint16 rl = (r[A] & 0x0f) - (x & 0x0f) - r.f.c;
  r[A] = rh;
  r.f.z = (uint8)rh == 0;
  r.f.n = 1;
  r.f.h = rl > 0x0f;
  r.f.c = rh > 0xff;
}

template<unsigned x> inline void CPU::op_sbc_a_r() { opi_sbc_a(r[x]); }
inline void CPU::op_sbc_a_n() { opi_sbc_a(op_read(r[PC]++)); }
inline void CPU::op_sbc_a_hl() { opi_sbc_a(op_read(r[HL])); }

inline void CPU::opi_and_a(uint8 x) {
  r[A] &= x;
  r.f.z = r[A] == 0;
  r.f.n = 0;
  r.f.h = 1;
  r.f.c = 0;
}

template<unsigned x> inline void CPU::op_and_a_r() { opi_and_a(r[x]); }
inline void CPU::op_and_a_n() { opi_and_a(op_read(r[PC]++)); }
inline void CPU::op_and_a_hl() { opi_and_a(op_read(r[HL])); }

inline void CPU::opi_xor_a(uint8 x) {
  r[A] ^= x;
  r.f.z = r[A] == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = 0;
}

template<unsigned x> inline void CPU::op_xor_a_r() { opi_xor_a(r[x]); }
inline void CPU::op_xor_a_n() { opi_xor_a(op_read(r[PC]++)); }
inline void CPU::op_xor_a_hl() { opi_xor_a(op_read(r[HL])); }

inline void CPU::opi_or_a(uint8 x) {
  r[A] |= x;
  r.f.z = r[A] == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = 0;
}

template<unsigned x> inline void CPU::op_or_a_r() { opi_or_a(r[x]); }
inline void CPU::op_or_a_n() { opi_or_a(op_read(r[PC]++)); }
inline void CPU::op_or_a_hl() { opi_or_a(op_read(r[HL])); }

inline void CPU::opi_cp_a(uint8 x) {
  uint16 rh = r[A] - x;
  uint16 rl = (r[A] & 0x0f) - (x & 0x0f);
  r.f.z = (uint8)rh == 0;
  r.f.n = 1;
  r.f.h = rl > 0x0f;
  r.f.c = rh > 0xff;
}

template<unsigned x> inline void CPU::op_cp_a_r() { opi_cp_a(r[x]); }
inline void CPU::op_cp_a_n() { opi_cp_a(op_read(r[PC]++)); }
inline void CPU::op_cp_a_hl() { opi_cp_a(op_read(r[HL])); }

template<unsigned x> inline void CPU::op_inc_r() {
  r[x]++;
  r.f.z = r[x] == 0;
  r.f.n = 0;
  r.f.h = (r[x] & 0x0f) == 0x00;
}

inline void CPU::op_inc_hl() {
  uint8 n = op_read(r[HL]);
  op_write(r[HL], ++n);
  r.f.z = n == 0;
  r.f.n = 0;
  r.f.h = (n & 0x0f) == 0x00;
}

template<unsigned x> inline void CPU::op_dec_r() {
  r[x]--;
  r.f.z = r[x] == 0;
  r.f.n = 1;
  r.f.h = (r[x] & 0x0f) == 0x0f;
}

inline void CPU::op_dec_hl() {
  uint8 n = op_read(r[HL]);
  op_write(r[HL], --n);
  r.f.z = n == 0;
  r.f.n = 1;
  r.f.h = (n & 0x0f) == 0x0f;
}

inline void CPU::op_daa() {
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

inline void CPU::op_cpl() {
  r[A] ^= 0xff;
  r.f.n = 1;
  r.f.h = 1;
}

//16-bit arithmetic commands

template<unsigned x> inline void CPU::op_add_hl_rr() {
  op_io();
  uint32 rb = (r[HL] + r[x]);
  uint32 rn = (r[HL] & 0xfff) + (r[x] & 0xfff);
  r[HL] = rb;
  r.f.n = 0;
  r.f.h = rn > 0x0fff;
  r.f.c = rb > 0xffff;
}

template<unsigned x> inline void CPU::op_inc_rr() {
  op_io();
  r[x]++;
}

template<unsigned x> inline void CPU::op_dec_rr() {
  op_io();
  r[x]--;
}

inline void CPU::op_add_sp_n() {
  op_io();
  op_io();
  signed n = (int8)op_read(r[PC]++);
  r.f.z = 0;
  r.f.n = 0;
  r.f.h = ((r[SP] & 0x0f) + (n & 0x0f)) > 0x0f;
  r.f.c = ((r[SP] & 0xff) + (n & 0xff)) > 0xff;
  r[SP] += n;
}

inline void CPU::op_ld_hl_sp_n() {
  op_io();
  signed n = (int8)op_read(r[PC]++);
  r.f.z = 0;
  r.f.n = 0;
  r.f.h = ((r[SP] & 0x0f) + (n & 0x0f)) > 0x0f;
  r.f.c = ((r[SP] & 0xff) + (n & 0xff)) > 0xff;
  r[HL] = r[SP] + n;
}

//rotate/shift commands

inline void CPU::op_rlca() {
  r[A] = (r[A] << 1) | (r[A] >> 7);
  r.f.z = 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = r[A] & 0x01;
}

inline void CPU::op_rla() {
  bool c = r[A] & 0x80;
  r[A] = (r[A] << 1) | (r.f.c << 0);
  r.f.z = 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = c;
}

inline void CPU::op_rrca() {
  r[A] = (r[A] >> 1) | (r[A] << 7);
  r.f.z = 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = r[A] & 0x80;
}

inline void CPU::op_rra() {
  bool c = r[A] & 0x01;
  r[A] = (r[A] >> 1) | (r.f.c << 7);
  r.f.z = 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = c;
}

template<unsigned x> inline void CPU::op_rlc_r() {
  r[x] = (r[x] << 1) | (r[x] >> 7);
  r.f.z = r[x] == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = r[x] & 0x01;
}

inline void CPU::op_rlc_hl() {
  uint8 n = op_read(r[HL]);
  n = (n << 1) | (n >> 7);
  op_write(r[HL], n);
  r.f.z = n == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = n & 0x01;
}

template<unsigned x> inline void CPU::op_rl_r() {
  bool c = r[x] & 0x80;
  r[x] = (r[x] << 1) | (r.f.c << 0);
  r.f.z = r[x] == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = c;
}

inline void CPU::op_rl_hl() {
  uint8 n = op_read(r[HL]);
  bool c = n & 0x80;
  n = (n << 1) | (r.f.c << 0);
  op_write(r[HL], n);
  r.f.z = n == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = c;
}

template<unsigned x> inline void CPU::op_rrc_r() {
  r[x] = (r[x] >> 1) | (r[x] << 7);
  r.f.z = r[x] == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = r[x] & 0x80;
}

inline void CPU::op_rrc_hl() {
  uint8 n = op_read(r[HL]);
  n = (n >> 1) | (n << 7);
  op_write(r[HL], n);
  r.f.z = n == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = n & 0x80;
}

template<unsigned x> inline void CPU::op_rr_r() {
  bool c = r[x] & 0x01;
  r[x] = (r[x] >> 1) | (r.f.c << 7);
  r.f.z = r[x] == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = c;
}

inline void CPU::op_rr_hl() {
  uint8 n = op_read(r[HL]);
  bool c = n & 0x01;
  n = (n >> 1) | (r.f.c << 7);
  op_write(r[HL], n);
  r.f.z = n == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = c;
}

template<unsigned x> inline void CPU::op_sla_r() {
  bool c = r[x] & 0x80;
  r[x] <<= 1;
  r.f.z = r[x] == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = c;
}

inline void CPU::op_sla_hl() {
  uint8 n = op_read(r[HL]);
  bool c = n & 0x80;
  n <<= 1;
  op_write(r[HL], n);
  r.f.z = n == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = c;
}

template<unsigned x> inline void CPU::op_swap_r() {
  r[x] = (r[x] << 4) | (r[x] >> 4);
  r.f.z = r[x] == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = 0;
}

inline void CPU::op_swap_hl() {
  uint8 n = op_read(r[HL]);
  n = (n << 4) | (n >> 4);
  op_write(r[HL], n);
  r.f.z = n == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = 0;
}

template<unsigned x> inline void CPU::op_sra_r() {
  bool c = r[x] & 0x01;
  r[x] = (int8)r[x] >> 1;
  r.f.z = r[x] == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = c;
}

inline void CPU::op_sra_hl() {
  uint8 n = op_read(r[HL]);
  bool c = n & 0x01;
  n = (int8)n >> 1;
  op_write(r[HL], n);
  r.f.z = n == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = c;
}

template<unsigned x> inline void CPU::op_srl_r() {
  bool c = r[x] & 0x01;
  r[x] >>= 1;
  r.f.z = r[x] == 0;
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = c;
}

inline void CPU::op_srl_hl() {
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

template<unsigned b, unsigned x> inline void CPU::op_bit_n_r() {
  r.f.z = (r[x] & (1 << b)) == 0;
  r.f.n = 0;
  r.f.h = 1;
}

template<unsigned b> inline void CPU::op_bit_n_hl() {
  uint8 n = op_read(r[HL]);
  r.f.z = (n & (1 << b)) == 0;
  r.f.n = 0;
  r.f.h = 1;
}

template<unsigned b, unsigned x> inline void CPU::op_set_n_r() {
  r[x] |= 1 << b;
}

template<unsigned b> inline void CPU::op_set_n_hl() {
  uint8 n = op_read(r[HL]);
  n |= 1 << b;
  op_write(r[HL], n);
}

template<unsigned b, unsigned x> inline void CPU::op_res_n_r() {
  r[x] &= ~(1 << b);
}

template<unsigned b> inline void CPU::op_res_n_hl() {
  uint8 n = op_read(r[HL]);
  n &= ~(1 << b);
  op_write(r[HL], n);
}

//control commands

inline void CPU::op_ccf() {
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = !r.f.c;
}

inline void CPU::op_scf() {
  r.f.n = 0;
  r.f.h = 0;
  r.f.c = 1;
}

inline void CPU::op_nop() {
}

inline void CPU::op_halt() {
  r.halt = true;
  while(r.halt == true) op_io();
}

inline void CPU::op_stop() {
  if(stop()) return;
  r.stop = true;
  while(r.stop == true) op_io();
}

inline void CPU::op_di() {
  r.ime = 0;
}

inline void CPU::op_ei() {
  r.ei = true;
//r.ime = 1;
}

//jump commands

inline void CPU::op_jp_nn() {
  uint8 lo = op_read(r[PC]++);
  uint8 hi = op_read(r[PC]++);
  r[PC] = (hi << 8) | (lo << 0);
  op_io();
}

inline void CPU::op_jp_hl() {
  r[PC] = r[HL];
}

template<unsigned x, bool y> inline void CPU::op_jp_f_nn() {
  uint8 lo = op_read(r[PC]++);
  uint8 hi = op_read(r[PC]++);
  if(r.f[x] == y) {
    r[PC] = (hi << 8) | (lo << 0);
    op_io();
  }
}

inline void CPU::op_jr_n() {
  int8 n = op_read(r[PC]++);
  r[PC] += n;
  op_io();
}

template<unsigned x, bool y> inline void CPU::op_jr_f_n() {
  int8 n = op_read(r[PC]++);
  if(r.f[x] == y) {
    r[PC] += n;
    op_io();
  }
}

inline void CPU::op_call_nn() {
  uint8 lo = op_read(r[PC]++);
  uint8 hi = op_read(r[PC]++);
  op_write(--r[SP], r[PC] >> 8);
  op_write(--r[SP], r[PC] >> 0);
  r[PC] = (hi << 8) | (lo << 0);
  op_io();
}

template<unsigned x, bool y> inline void CPU::op_call_f_nn() {
  uint8 lo = op_read(r[PC]++);
  uint8 hi = op_read(r[PC]++);
  if(r.f[x] == y) {
    op_write(--r[SP], r[PC] >> 8);
    op_write(--r[SP], r[PC] >> 0);
    r[PC] = (hi << 8) | (lo << 0);
    op_io();
  }
}

inline void CPU::op_ret() {
  uint8 lo = op_read(r[SP]++);
  uint8 hi = op_read(r[SP]++);
  r[PC] = (hi << 8) | (lo << 0);
  op_io();
}

template<unsigned x, bool y> inline void CPU::op_ret_f() {
  op_io();
  if(r.f[x] == y) {
    uint8 lo = op_read(r[SP]++);
    uint8 hi = op_read(r[SP]++);
    r[PC] = (hi << 8) | (lo << 0);
    op_io();
  }
}

inline void CPU::op_reti() {
  uint8 lo = op_read(r[SP]++);
  uint8 hi = op_read(r[SP]++);
  r[PC] = (hi << 8) | (lo << 0);
  op_io();
  r.ime = 1;
}

template<unsigned n> inline void CPU::op_rst_n() {
  op_write(--r[SP], r[PC] >> 8);
  op_write(--r[SP], r[PC] >> 0);
  r[PC] = n;
  op_io();
}

void CPU::dump() {
    std::ofstream txt, csv;
    txt.open("/tmp/benchmarks.txt", std::ios::out);

    txt  << "Benchmarks for Termboy:" << std::endl;
    txt  << "-------------------------------------------" << std::endl;
    txt  << "\tTotal Instructions:\t\t" << instruction_count << std::endl;
    txt  << "\tTotal Time:\t\t\t" << ((double) std::clock() - global_time ) / 1000 << std::endl;
    txt  << "\tMax Time:\t\t\t" << max_time << std::endl;
    txt  << "\tTime/Instruction:\t\t" << time << std::endl;
    txt  << "\tAvg Time in synch:\t\t\t" << synch_time << std::endl;

    txt.close();
}

void CPU::get_times(double final_time, double synch) {
    if(max_time < final_time )
        max_time = final_time;

    time = (final_time + time) / 2;
    synch_time = (synch_time + synch)/2;
}



void CPU::exec_cb() { 
#include "labels.cpp"

 uint8 opcode = op_read(r[PC]++);  
 cb_operation=true;
 goto *labels[opcode];

l00: CPU::op_rlc_r<B>(); goto end;
l01: CPU::op_rlc_r<C>(); goto end;
l02: CPU::op_rlc_r<D>(); goto end;
l03: CPU::op_rlc_r<E>(); goto end;
l04: CPU::op_rlc_r<H>(); goto end;
l05: CPU::op_rlc_r<L>(); goto end;
l06: CPU::op_rlc_hl(); goto end;
l07: CPU::op_rlc_r<A>(); goto end;
l08: CPU::op_rrc_r<B>(); goto end;
l09: CPU::op_rrc_r<C>(); goto end;
l0a: CPU::op_rrc_r<D>(); goto end;
l0b: CPU::op_rrc_r<E>(); goto end;
l0c: CPU::op_rrc_r<H>(); goto end;
l0d: CPU::op_rrc_r<L>(); goto end;
l0e: CPU::op_rrc_hl(); goto end;
l0f: CPU::op_rrc_r<A>(); goto end;
l10: CPU::op_rl_r<B>(); goto end;
l11: CPU::op_rl_r<C>(); goto end;
l12: CPU::op_rl_r<D>(); goto end;
l13: CPU::op_rl_r<E>(); goto end;
l14: CPU::op_rl_r<H>(); goto end;
l15: CPU::op_rl_r<L>(); goto end;
l16: CPU::op_rl_hl(); goto end;
l17: CPU::op_rl_r<A>(); goto end;
l18: CPU::op_rr_r<B>(); goto end;
l19: CPU::op_rr_r<C>(); goto end;
l1a: CPU::op_rr_r<D>(); goto end;
l1b: CPU::op_rr_r<E>(); goto end;
l1c: CPU::op_rr_r<H>(); goto end;
l1d: CPU::op_rr_r<L>(); goto end;
l1e: CPU::op_rr_hl(); goto end;
l1f: CPU::op_rr_r<A>(); goto end;
l20: CPU::op_sla_r<B>(); goto end;
l21: CPU::op_sla_r<C>(); goto end;
l22: CPU::op_sla_r<D>(); goto end;
l23: CPU::op_sla_r<E>(); goto end;
l24: CPU::op_sla_r<H>(); goto end;
l25: CPU::op_sla_r<L>(); goto end;
l26: CPU::op_sla_hl(); goto end;
l27: CPU::op_sla_r<A>(); goto end;
l28: CPU::op_sra_r<B>(); goto end;
l29: CPU::op_sra_r<C>(); goto end;
l2a: CPU::op_sra_r<D>(); goto end;
l2b: CPU::op_sra_r<E>(); goto end;
l2c: CPU::op_sra_r<H>(); goto end;
l2d: CPU::op_sra_r<L>(); goto end;
l2e: CPU::op_sra_hl(); goto end;
l2f: CPU::op_sra_r<A>(); goto end;
l30: CPU::op_swap_r<B>(); goto end;
l31: CPU::op_swap_r<C>(); goto end;
l32: CPU::op_swap_r<D>(); goto end;
l33: CPU::op_swap_r<E>(); goto end;
l34: CPU::op_swap_r<H>(); goto end;
l35: CPU::op_swap_r<L>(); goto end;
l36: CPU::op_swap_hl(); goto end;
l37: CPU::op_swap_r<A>(); goto end;
l38: CPU::op_srl_r<B>(); goto end;
l39: CPU::op_srl_r<C>(); goto end;
l3a: CPU::op_srl_r<D>(); goto end;
l3b: CPU::op_srl_r<E>(); goto end;
l3c: CPU::op_srl_r<H>(); goto end;
l3d: CPU::op_srl_r<L>(); goto end;
l3e: CPU::op_srl_hl(); goto end;
l3f: CPU::op_srl_r<A>(); goto end;
l40: CPU::op_bit_n_r<0, B>(); goto end;
l41: CPU::op_bit_n_r<0, C>(); goto end;
l42: CPU::op_bit_n_r<0, D>(); goto end;
l43: CPU::op_bit_n_r<0, E>(); goto end;
l44: CPU::op_bit_n_r<0, H>(); goto end;
l45: CPU::op_bit_n_r<0, L>(); goto end;
l46: CPU::op_bit_n_hl<0>(); goto end;
l47: CPU::op_bit_n_r<0, A>(); goto end;
l48: CPU::op_bit_n_r<1, B>(); goto end;
l49: CPU::op_bit_n_r<1, C>(); goto end;
l4a: CPU::op_bit_n_r<1, D>(); goto end;
l4b: CPU::op_bit_n_r<1, E>(); goto end;
l4c: CPU::op_bit_n_r<1, H>(); goto end;
l4d: CPU::op_bit_n_r<1, L>(); goto end;
l4e: CPU::op_bit_n_hl<1>(); goto end;
l4f: CPU::op_bit_n_r<1, A>(); goto end;
l50: CPU::op_bit_n_r<2, B>(); goto end;
l51: CPU::op_bit_n_r<2, C>(); goto end;
l52: CPU::op_bit_n_r<2, D>(); goto end;
l53: CPU::op_bit_n_r<2, E>(); goto end;
l54: CPU::op_bit_n_r<2, H>(); goto end;
l55: CPU::op_bit_n_r<2, L>(); goto end;
l56: CPU::op_bit_n_hl<2>(); goto end;
l57: CPU::op_bit_n_r<2, A>(); goto end;
l58: CPU::op_bit_n_r<3, B>(); goto end;
l59: CPU::op_bit_n_r<3, C>(); goto end;
l5a: CPU::op_bit_n_r<3, D>(); goto end;
l5b: CPU::op_bit_n_r<3, E>(); goto end;
l5c: CPU::op_bit_n_r<3, H>(); goto end;
l5d: CPU::op_bit_n_r<3, L>(); goto end;
l5e: CPU::op_bit_n_hl<3>(); goto end;
l5f: CPU::op_bit_n_r<3, A>(); goto end;
l60: CPU::op_bit_n_r<4, B>(); goto end;
l61: CPU::op_bit_n_r<4, C>(); goto end;
l62: CPU::op_bit_n_r<4, D>(); goto end;
l63: CPU::op_bit_n_r<4, E>(); goto end;
l64: CPU::op_bit_n_r<4, H>(); goto end;
l65: CPU::op_bit_n_r<4, L>(); goto end;
l66: CPU::op_bit_n_hl<4>(); goto end;
l67: CPU::op_bit_n_r<4, A>(); goto end;
l68: CPU::op_bit_n_r<5, B>(); goto end;
l69: CPU::op_bit_n_r<5, C>(); goto end;
l6a: CPU::op_bit_n_r<5, D>(); goto end;
l6b: CPU::op_bit_n_r<5, E>(); goto end;
l6c: CPU::op_bit_n_r<5, H>(); goto end;
l6d: CPU::op_bit_n_r<5, L>(); goto end;
l6e: CPU::op_bit_n_hl<5>(); goto end;
l6f: CPU::op_bit_n_r<5, A>(); goto end;
l70: CPU::op_bit_n_r<6, B>(); goto end;
l71: CPU::op_bit_n_r<6, C>(); goto end;
l72: CPU::op_bit_n_r<6, D>(); goto end;
l73: CPU::op_bit_n_r<6, E>(); goto end;
l74: CPU::op_bit_n_r<6, H>(); goto end;
l75: CPU::op_bit_n_r<6, L>(); goto end;
l76: CPU::op_bit_n_hl<6>(); goto end;
l77: CPU::op_bit_n_r<6, A>(); goto end;
l78: CPU::op_bit_n_r<7, B>(); goto end;
l79: CPU::op_bit_n_r<7, C>(); goto end;
l7a: CPU::op_bit_n_r<7, D>(); goto end;
l7b: CPU::op_bit_n_r<7, E>(); goto end;
l7c: CPU::op_bit_n_r<7, H>(); goto end;
l7d: CPU::op_bit_n_r<7, L>(); goto end;
l7e: CPU::op_bit_n_hl<7>(); goto end;
l7f: CPU::op_bit_n_r<7, A>(); goto end;
l80: CPU::op_res_n_r<0, B>(); goto end;
l81: CPU::op_res_n_r<0, C>(); goto end;
l82: CPU::op_res_n_r<0, D>(); goto end;
l83: CPU::op_res_n_r<0, E>(); goto end;
l84: CPU::op_res_n_r<0, H>(); goto end;
l85: CPU::op_res_n_r<0, L>(); goto end;
l86: CPU::op_res_n_hl<0>(); goto end;
l87: CPU::op_res_n_r<0, A>(); goto end;
l88: CPU::op_res_n_r<1, B>(); goto end;
l89: CPU::op_res_n_r<1, C>(); goto end;
l8a: CPU::op_res_n_r<1, D>(); goto end;
l8b: CPU::op_res_n_r<1, E>(); goto end;
l8c: CPU::op_res_n_r<1, H>(); goto end;
l8d: CPU::op_res_n_r<1, L>(); goto end;
l8e: CPU::op_res_n_hl<1>(); goto end;
l8f: CPU::op_res_n_r<1, A>(); goto end;
l90: CPU::op_res_n_r<2, B>(); goto end;
l91: CPU::op_res_n_r<2, C>(); goto end;
l92: CPU::op_res_n_r<2, D>(); goto end;
l93: CPU::op_res_n_r<2, E>(); goto end;
l94: CPU::op_res_n_r<2, H>(); goto end;
l95: CPU::op_res_n_r<2, L>(); goto end;
l96: CPU::op_res_n_hl<2>(); goto end;
l97: CPU::op_res_n_r<2, A>(); goto end;
l98: CPU::op_res_n_r<3, B>(); goto end;
l99: CPU::op_res_n_r<3, C>(); goto end;
l9a: CPU::op_res_n_r<3, D>(); goto end;
l9b: CPU::op_res_n_r<3, E>(); goto end;
l9c: CPU::op_res_n_r<3, H>(); goto end;
l9d: CPU::op_res_n_r<3, L>(); goto end;
l9e: CPU::op_res_n_hl<3>(); goto end;
l9f: CPU::op_res_n_r<3, A>(); goto end;
la0: CPU::op_res_n_r<4, B>(); goto end;
la1: CPU::op_res_n_r<4, C>(); goto end;
la2: CPU::op_res_n_r<4, D>(); goto end;
la3: CPU::op_res_n_r<4, E>(); goto end;
la4: CPU::op_res_n_r<4, H>(); goto end;
la5: CPU::op_res_n_r<4, L>(); goto end;
la6: CPU::op_res_n_hl<4>(); goto end;
la7: CPU::op_res_n_r<4, A>(); goto end;
la8: CPU::op_res_n_r<5, B>(); goto end;
la9: CPU::op_res_n_r<5, C>(); goto end;
laa: CPU::op_res_n_r<5, D>(); goto end;
lab: CPU::op_res_n_r<5, E>(); goto end;
lac: CPU::op_res_n_r<5, H>(); goto end;
lad: CPU::op_res_n_r<5, L>(); goto end;
lae: CPU::op_res_n_hl<5>(); goto end;
laf: CPU::op_res_n_r<5, A>(); goto end;
lb0: CPU::op_res_n_r<6, B>(); goto end;
lb1: CPU::op_res_n_r<6, C>(); goto end;
lb2: CPU::op_res_n_r<6, D>(); goto end;
lb3: CPU::op_res_n_r<6, E>(); goto end;
lb4: CPU::op_res_n_r<6, H>(); goto end;
lb5: CPU::op_res_n_r<6, L>(); goto end;
lb6: CPU::op_res_n_hl<6>(); goto end;
lb7: CPU::op_res_n_r<6, A>(); goto end;
lb8: CPU::op_res_n_r<7, B>(); goto end;
lb9: CPU::op_res_n_r<7, C>(); goto end;
lba: CPU::op_res_n_r<7, D>(); goto end;
lbb: CPU::op_res_n_r<7, E>(); goto end;
lbc: CPU::op_res_n_r<7, H>(); goto end;
lbd: CPU::op_res_n_r<7, L>(); goto end;
lbe: CPU::op_res_n_hl<7>(); goto end;
lbf: CPU::op_res_n_r<7, A>(); goto end;
lc0: CPU::op_set_n_r<0, B>(); goto end;
lc1: CPU::op_set_n_r<0, C>(); goto end;
lc2: CPU::op_set_n_r<0, D>(); goto end;
lc3: CPU::op_set_n_r<0, E>(); goto end;
lc4: CPU::op_set_n_r<0, H>(); goto end;
lc5: CPU::op_set_n_r<0, L>(); goto end;
lc6: CPU::op_set_n_hl<0>(); goto end;
lc7: CPU::op_set_n_r<0, A>(); goto end;
lc8: CPU::op_set_n_r<1, B>(); goto end;
lc9: CPU::op_set_n_r<1, C>(); goto end;
lca: CPU::op_set_n_r<1, D>(); goto end;
lcb: CPU::op_set_n_r<1, E>(); goto end;
lcc: CPU::op_set_n_r<1, H>(); goto end;
lcd: CPU::op_set_n_r<1, L>(); goto end;
lce: CPU::op_set_n_hl<1>(); goto end;
lcf: CPU::op_set_n_r<1, A>(); goto end;
ld0: CPU::op_set_n_r<2, B>(); goto end;
ld1: CPU::op_set_n_r<2, C>(); goto end;
ld2: CPU::op_set_n_r<2, D>(); goto end;
ld3: CPU::op_set_n_r<2, E>(); goto end;
ld4: CPU::op_set_n_r<2, H>(); goto end;
ld5: CPU::op_set_n_r<2, L>(); goto end;
ld6: CPU::op_set_n_hl<2>(); goto end;
ld7: CPU::op_set_n_r<2, A>(); goto end;
ld8: CPU::op_set_n_r<3, B>(); goto end;
ld9: CPU::op_set_n_r<3, C>(); goto end;
lda: CPU::op_set_n_r<3, D>(); goto end;
ldb: CPU::op_set_n_r<3, E>(); goto end;
ldc: CPU::op_set_n_r<3, H>(); goto end;
ldd: CPU::op_set_n_r<3, L>(); goto end;
lde: CPU::op_set_n_hl<3>(); goto end;
ldf: CPU::op_set_n_r<3, A>(); goto end;
le0: CPU::op_set_n_r<4, B>(); goto end;
le1: CPU::op_set_n_r<4, C>(); goto end;
le2: CPU::op_set_n_r<4, D>(); goto end;
le3: CPU::op_set_n_r<4, E>(); goto end;
le4: CPU::op_set_n_r<4, H>(); goto end;
le5: CPU::op_set_n_r<4, L>(); goto end;
le6: CPU::op_set_n_hl<4>(); goto end;
le7: CPU::op_set_n_r<4, A>(); goto end;
le8: CPU::op_set_n_r<5, B>(); goto end;
le9: CPU::op_set_n_r<5, C>(); goto end;
lea: CPU::op_set_n_r<5, D>(); goto end;
leb: CPU::op_set_n_r<5, E>(); goto end;
lec: CPU::op_set_n_r<5, H>(); goto end;
led: CPU::op_set_n_r<5, L>(); goto end;
lee: CPU::op_set_n_hl<5>(); goto end;
lef: CPU::op_set_n_r<5, A>(); goto end;
lf0: CPU::op_set_n_r<6, B>(); goto end;
lf1: CPU::op_set_n_r<6, C>(); goto end;
lf2: CPU::op_set_n_r<6, D>(); goto end;
lf3: CPU::op_set_n_r<6, E>(); goto end;
lf4: CPU::op_set_n_r<6, H>(); goto end;
lf5: CPU::op_set_n_r<6, L>(); goto end;
lf6: CPU::op_set_n_hl<6>(); goto end;
lf7: CPU::op_set_n_r<6, A>(); goto end;
lf8: CPU::op_set_n_r<7, B>(); goto end;
lf9: CPU::op_set_n_r<7, C>(); goto end;
lfa: CPU::op_set_n_r<7, D>(); goto end;
lfb: CPU::op_set_n_r<7, E>(); goto end;
lfc: CPU::op_set_n_r<7, H>(); goto end;
lfd: CPU::op_set_n_r<7, L>(); goto end;
lfe: CPU::op_set_n_hl<7>(); goto end;
lff: CPU::op_set_n_r<7, A>(); goto end;

end: return;
}

#define d_dispatcher() \
  final_time = double(std::clock() - exec_time) / (double)CLOCKS_PER_SEC; \
  get_times(final_time, synch_time); ddispatcher()

#define ddispatcher() \
  synch_time = std::clock(); \
  if(scheduler.sync == Scheduler::SynchronizeMode::CPU) { \
     scheduler.sync = Scheduler::SynchronizeMode::All;    \
     scheduler.exit(Scheduler::ExitReason::SynchronizeEvent); \
  } test_for_interrupt(); \
  synch_time = double(std::clock() - synch_time) / (double)CLOCKS_PER_SEC; \
  exec_time = std::clock(); \
  opcode = op_read(r[PC]++); \
  instruction_count++; \
  goto *labels[opcode]

void CPU::exec() {
#include "labels.cpp"

  std::clock_t exec_time;
  double final_time;
  uint8 opcode; 

  global_time = std::clock();
  ddispatcher();


l00: CPU::op_nop(); d_dispatcher();
l01: CPU::op_ld_rr_nn<BC>(); d_dispatcher();
l02: CPU::op_ld_rr_a<BC>(); d_dispatcher();
l03: CPU::op_inc_rr<BC>(); d_dispatcher();
l04: CPU::op_inc_r<B>(); d_dispatcher();
l05: CPU::op_dec_r<B>(); d_dispatcher();
l06: CPU::op_ld_r_n<B>(); d_dispatcher();
l07: CPU::op_rlca(); d_dispatcher();
l08: CPU::op_ld_nn_sp(); d_dispatcher();
l09: CPU::op_add_hl_rr<BC>(); d_dispatcher();
l0a: CPU::op_ld_a_rr<BC>(); d_dispatcher();
l0b: CPU::op_dec_rr<BC>(); d_dispatcher();
l0c: CPU::op_inc_r<C>(); d_dispatcher();
l0d: CPU::op_dec_r<C>(); d_dispatcher();
l0e: CPU::op_ld_r_n<C>(); d_dispatcher();
l0f: CPU::op_rrca(); d_dispatcher();
l10: CPU::op_stop(); d_dispatcher();
l11: CPU::op_ld_rr_nn<DE>(); d_dispatcher();
l12: CPU::op_ld_rr_a<DE>(); d_dispatcher();
l13: CPU::op_inc_rr<DE>(); d_dispatcher();
l14: CPU::op_inc_r<D>(); d_dispatcher();
l15: CPU::op_dec_r<D>(); d_dispatcher();
l16: CPU::op_ld_r_n<D>(); d_dispatcher();
l17: CPU::op_rla(); d_dispatcher();
l18: CPU::op_jr_n(); d_dispatcher();
l19: CPU::op_add_hl_rr<DE>(); d_dispatcher();
l1a: CPU::op_ld_a_rr<DE>(); d_dispatcher();
l1b: CPU::op_dec_rr<DE>(); d_dispatcher();
l1c: CPU::op_inc_r<E>(); d_dispatcher();
l1d: CPU::op_dec_r<E>(); d_dispatcher();
l1e: CPU::op_ld_r_n<E>(); d_dispatcher();
l1f: CPU::op_rra(); d_dispatcher();
l20: CPU::op_jr_f_n<ZF, 0>(); d_dispatcher();
l21: CPU::op_ld_rr_nn<HL>(); d_dispatcher();
l22: CPU::op_ldi_hl_a(); d_dispatcher();
l23: CPU::op_inc_rr<HL>(); d_dispatcher();
l24: CPU::op_inc_r<H>(); d_dispatcher();
l25: CPU::op_dec_r<H>(); d_dispatcher();
l26: CPU::op_ld_r_n<H>(); d_dispatcher();
l27: CPU::op_daa(); d_dispatcher();
l28: CPU::op_jr_f_n<ZF, 1>(); d_dispatcher();
l29: CPU::op_add_hl_rr<HL>(); d_dispatcher();
l2a: CPU::op_ldi_a_hl(); d_dispatcher();
l2b: CPU::op_dec_rr<HL>(); d_dispatcher();
l2c: CPU::op_inc_r<L>(); d_dispatcher();
l2d: CPU::op_dec_r<L>(); d_dispatcher();
l2e: CPU::op_ld_r_n<L>(); d_dispatcher();
l2f: CPU::op_cpl(); d_dispatcher();
l30: CPU::op_jr_f_n<CF, 0>(); d_dispatcher();
l31: CPU::op_ld_rr_nn<SP>(); d_dispatcher();
l32: CPU::op_ldd_hl_a(); d_dispatcher();
l33: CPU::op_inc_rr<SP>(); d_dispatcher();
l34: CPU::op_inc_hl(); d_dispatcher();
l35: CPU::op_dec_hl(); d_dispatcher();
l36: CPU::op_ld_hl_n(); d_dispatcher();
l37: CPU::op_scf(); d_dispatcher();
l38: CPU::op_jr_f_n<CF, 1>(); d_dispatcher();
l39: CPU::op_add_hl_rr<SP>(); d_dispatcher();
l3a: CPU::op_ldd_a_hl(); d_dispatcher();
l3b: CPU::op_dec_rr<SP>(); d_dispatcher();
l3c: CPU::op_inc_r<A>(); d_dispatcher();
l3d: CPU::op_dec_r<A>(); d_dispatcher();
l3e: CPU::op_ld_r_n<A>(); d_dispatcher();
l3f: CPU::op_ccf(); d_dispatcher();
l40: CPU::op_ld_r_r<B, B>(); d_dispatcher();
l41: CPU::op_ld_r_r<B, C>(); d_dispatcher();
l42: CPU::op_ld_r_r<B, D>(); d_dispatcher();
l43: CPU::op_ld_r_r<B, E>(); d_dispatcher();
l44: CPU::op_ld_r_r<B, H>(); d_dispatcher();
l45: CPU::op_ld_r_r<B, L>(); d_dispatcher();
l46: CPU::op_ld_r_hl<B>(); d_dispatcher();
l47: CPU::op_ld_r_r<B, A>(); d_dispatcher();
l48: CPU::op_ld_r_r<C, B>(); d_dispatcher();
l49: CPU::op_ld_r_r<C, C>(); d_dispatcher();
l4a: CPU::op_ld_r_r<C, D>(); d_dispatcher();
l4b: CPU::op_ld_r_r<C, E>(); d_dispatcher();
l4c: CPU::op_ld_r_r<C, H>(); d_dispatcher();
l4d: CPU::op_ld_r_r<C, L>(); d_dispatcher();
l4e: CPU::op_ld_r_hl<C>(); d_dispatcher();
l4f: CPU::op_ld_r_r<C, A>(); d_dispatcher();
l50: CPU::op_ld_r_r<D, B>(); d_dispatcher();
l51: CPU::op_ld_r_r<D, C>(); d_dispatcher();
l52: CPU::op_ld_r_r<D, D>(); d_dispatcher();
l53: CPU::op_ld_r_r<D, E>(); d_dispatcher();
l54: CPU::op_ld_r_r<D, H>(); d_dispatcher();
l55: CPU::op_ld_r_r<D, L>(); d_dispatcher();
l56: CPU::op_ld_r_hl<D>(); d_dispatcher();
l57: CPU::op_ld_r_r<D, A>(); d_dispatcher();
l58: CPU::op_ld_r_r<E, B>(); d_dispatcher();
l59: CPU::op_ld_r_r<E, C>(); d_dispatcher();
l5a: CPU::op_ld_r_r<E, D>(); d_dispatcher();
l5b: CPU::op_ld_r_r<E, E>(); d_dispatcher();
l5c: CPU::op_ld_r_r<E, H>(); d_dispatcher();
l5d: CPU::op_ld_r_r<E, L>(); d_dispatcher();
l5e: CPU::op_ld_r_hl<E>(); d_dispatcher();
l5f: CPU::op_ld_r_r<E, A>(); d_dispatcher();
l60: CPU::op_ld_r_r<H, B>(); d_dispatcher();
l61: CPU::op_ld_r_r<H, C>(); d_dispatcher();
l62: CPU::op_ld_r_r<H, D>(); d_dispatcher();
l63: CPU::op_ld_r_r<H, E>(); d_dispatcher();
l64: CPU::op_ld_r_r<H, H>(); d_dispatcher();
l65: CPU::op_ld_r_r<H, L>(); d_dispatcher();
l66: CPU::op_ld_r_hl<H>(); d_dispatcher();
l67: CPU::op_ld_r_r<H, A>(); d_dispatcher();
l68: CPU::op_ld_r_r<L, B>(); d_dispatcher();
l69: CPU::op_ld_r_r<L, C>(); d_dispatcher();
l6a: CPU::op_ld_r_r<L, D>(); d_dispatcher();
l6b: CPU::op_ld_r_r<L, E>(); d_dispatcher();
l6c: CPU::op_ld_r_r<L, H>(); d_dispatcher();
l6d: CPU::op_ld_r_r<L, L>(); d_dispatcher();
l6e: CPU::op_ld_r_hl<L>(); d_dispatcher();
l6f: CPU::op_ld_r_r<L, A>(); d_dispatcher();
l70: CPU::op_ld_hl_r<B>(); d_dispatcher();
l71: CPU::op_ld_hl_r<C>(); d_dispatcher();
l72: CPU::op_ld_hl_r<D>(); d_dispatcher();
l73: CPU::op_ld_hl_r<E>(); d_dispatcher();
l74: CPU::op_ld_hl_r<H>(); d_dispatcher();
l75: CPU::op_ld_hl_r<L>(); d_dispatcher();
l76: CPU::op_halt(); d_dispatcher();
l77: CPU::op_ld_hl_r<A>(); d_dispatcher();
l78: CPU::op_ld_r_r<A, B>(); d_dispatcher();
l79: CPU::op_ld_r_r<A, C>(); d_dispatcher();
l7a: CPU::op_ld_r_r<A, D>(); d_dispatcher();
l7b: CPU::op_ld_r_r<A, E>(); d_dispatcher();
l7c: CPU::op_ld_r_r<A, H>(); d_dispatcher();
l7d: CPU::op_ld_r_r<A, L>(); d_dispatcher();
l7e: CPU::op_ld_r_hl<A>(); d_dispatcher();
l7f: CPU::op_ld_r_r<A, A>(); d_dispatcher();
l80: CPU::op_add_a_r<B>(); d_dispatcher();
l81: CPU::op_add_a_r<C>(); d_dispatcher();
l82: CPU::op_add_a_r<D>(); d_dispatcher();
l83: CPU::op_add_a_r<E>(); d_dispatcher();
l84: CPU::op_add_a_r<H>(); d_dispatcher();
l85: CPU::op_add_a_r<L>(); d_dispatcher();
l86: CPU::op_add_a_hl(); d_dispatcher();
l87: CPU::op_add_a_r<A>(); d_dispatcher();
l88: CPU::op_adc_a_r<B>(); d_dispatcher();
l89: CPU::op_adc_a_r<C>(); d_dispatcher();
l8a: CPU::op_adc_a_r<D>(); d_dispatcher();
l8b: CPU::op_adc_a_r<E>(); d_dispatcher();
l8c: CPU::op_adc_a_r<H>(); d_dispatcher();
l8d: CPU::op_adc_a_r<L>(); d_dispatcher();
l8e: CPU::op_adc_a_hl(); d_dispatcher();
l8f: CPU::op_adc_a_r<A>(); d_dispatcher();
l90: CPU::op_sub_a_r<B>(); d_dispatcher();
l91: CPU::op_sub_a_r<C>(); d_dispatcher();
l92: CPU::op_sub_a_r<D>(); d_dispatcher();
l93: CPU::op_sub_a_r<E>(); d_dispatcher();
l94: CPU::op_sub_a_r<H>(); d_dispatcher();
l95: CPU::op_sub_a_r<L>(); d_dispatcher();
l96: CPU::op_sub_a_hl(); d_dispatcher();
l97: CPU::op_sub_a_r<A>(); d_dispatcher();
l98: CPU::op_sbc_a_r<B>(); d_dispatcher();
l99: CPU::op_sbc_a_r<C>(); d_dispatcher();
l9a: CPU::op_sbc_a_r<D>(); d_dispatcher();
l9b: CPU::op_sbc_a_r<E>(); d_dispatcher();
l9c: CPU::op_sbc_a_r<H>(); d_dispatcher();
l9d: CPU::op_sbc_a_r<L>(); d_dispatcher();
l9e: CPU::op_sbc_a_hl(); d_dispatcher();
l9f: CPU::op_sbc_a_r<A>(); d_dispatcher();
la0: CPU::op_and_a_r<B>(); d_dispatcher();
la1: CPU::op_and_a_r<C>(); d_dispatcher();
la2: CPU::op_and_a_r<D>(); d_dispatcher();
la3: CPU::op_and_a_r<E>(); d_dispatcher();
la4: CPU::op_and_a_r<H>(); d_dispatcher();
la5: CPU::op_and_a_r<L>(); d_dispatcher();
la6: CPU::op_and_a_hl(); d_dispatcher();
la7: CPU::op_and_a_r<A>(); d_dispatcher();
la8: CPU::op_xor_a_r<B>(); d_dispatcher();
la9: CPU::op_xor_a_r<C>(); d_dispatcher();
laa: CPU::op_xor_a_r<D>(); d_dispatcher();
lab: CPU::op_xor_a_r<E>(); d_dispatcher();
lac: CPU::op_xor_a_r<H>(); d_dispatcher();
lad: CPU::op_xor_a_r<L>(); d_dispatcher();
lae: CPU::op_xor_a_hl(); d_dispatcher();
laf: CPU::op_xor_a_r<A>(); d_dispatcher();
lb0: CPU::op_or_a_r<B>(); d_dispatcher();
lb1: CPU::op_or_a_r<C>(); d_dispatcher();
lb2: CPU::op_or_a_r<D>(); d_dispatcher();
lb3: CPU::op_or_a_r<E>(); d_dispatcher();
lb4: CPU::op_or_a_r<H>(); d_dispatcher();
lb5: CPU::op_or_a_r<L>(); d_dispatcher();
lb6: CPU::op_or_a_hl(); d_dispatcher();
lb7: CPU::op_or_a_r<A>(); d_dispatcher();
lb8: CPU::op_cp_a_r<B>(); d_dispatcher();
lb9: CPU::op_cp_a_r<C>(); d_dispatcher();
lba: CPU::op_cp_a_r<D>(); d_dispatcher();
lbb: CPU::op_cp_a_r<E>(); d_dispatcher();
lbc: CPU::op_cp_a_r<H>(); d_dispatcher();
lbd: CPU::op_cp_a_r<L>(); d_dispatcher();
lbe: CPU::op_cp_a_hl(); d_dispatcher();
lbf: CPU::op_cp_a_r<A>(); d_dispatcher();
lc0: CPU::op_ret_f<ZF, 0>(); d_dispatcher();
lc1: CPU::op_pop_rr<BC>(); d_dispatcher();
lc2: CPU::op_jp_f_nn<ZF, 0>(); d_dispatcher();
lc3: CPU::op_jp_nn(); d_dispatcher();
lc4: CPU::op_call_f_nn<ZF, 0>(); d_dispatcher();
lc5: CPU::op_push_rr<BC>(); d_dispatcher();
lc6: CPU::op_add_a_n(); d_dispatcher();
lc7: CPU::op_rst_n<0x00>(); d_dispatcher();
lc8: CPU::op_ret_f<ZF, 1>(); d_dispatcher();
lc9: CPU::op_ret(); d_dispatcher();
lca: CPU::op_jp_f_nn<ZF, 1>(); d_dispatcher();
lcb: CPU::op_cb(); d_dispatcher();
lcc: CPU::op_call_f_nn<ZF, 1>(); d_dispatcher();
lcd: CPU::op_call_nn(); d_dispatcher();
lce: CPU::op_adc_a_n(); d_dispatcher();
lcf: CPU::op_rst_n<0x08>(); d_dispatcher();
ld0: CPU::op_ret_f<CF, 0>(); d_dispatcher();
ld1: CPU::op_pop_rr<DE>(); d_dispatcher();
ld2: CPU::op_jp_f_nn<CF, 0>(); d_dispatcher();
ld3: CPU::op_xx(); d_dispatcher();
ld4: CPU::op_call_f_nn<CF, 0>(); d_dispatcher();
ld5: CPU::op_push_rr<DE>(); d_dispatcher();
ld6: CPU::op_sub_a_n(); d_dispatcher();
ld7: CPU::op_rst_n<0x10>(); d_dispatcher();
ld8: CPU::op_ret_f<CF, 1>(); d_dispatcher();
ld9: CPU::op_reti(); d_dispatcher();
lda: CPU::op_jp_f_nn<CF, 1>(); d_dispatcher();
ldb: CPU::op_xx(); d_dispatcher();
ldc: CPU::op_call_f_nn<CF, 1>(); d_dispatcher();
ldd: CPU::op_xx(); d_dispatcher();
lde: CPU::op_sbc_a_n(); d_dispatcher();
ldf: CPU::op_rst_n<0x18>(); d_dispatcher();
le0: CPU::op_ld_ffn_a(); d_dispatcher();
le1: CPU::op_pop_rr<HL>(); d_dispatcher();
le2: CPU::op_ld_ffc_a(); d_dispatcher();
le3: CPU::op_xx(); d_dispatcher();
le4: CPU::op_xx(); d_dispatcher();
le5: CPU::op_push_rr<HL>(); d_dispatcher();
le6: CPU::op_and_a_n(); d_dispatcher();
le7: CPU::op_rst_n<0x20>(); d_dispatcher();
le8: CPU::op_add_sp_n(); d_dispatcher();
le9: CPU::op_jp_hl(); d_dispatcher();
lea: CPU::op_ld_nn_a(); d_dispatcher();
leb: CPU::op_xx(); d_dispatcher();
lec: CPU::op_xx(); d_dispatcher();
led: CPU::op_xx(); d_dispatcher();
lee: CPU::op_xor_a_n(); d_dispatcher();
lef: CPU::op_rst_n<0x28>(); d_dispatcher();
lf0: CPU::op_ld_a_ffn(); d_dispatcher();
lf1: CPU::op_pop_rr<AF>(); d_dispatcher();
lf2: CPU::op_ld_a_ffc(); d_dispatcher();
lf3: CPU::op_di(); d_dispatcher();
lf4: CPU::op_xx(); d_dispatcher();
lf5: CPU::op_push_rr<AF>(); d_dispatcher();
lf6: CPU::op_or_a_n(); d_dispatcher();
lf7: CPU::op_rst_n<0x30>(); d_dispatcher();
lf8: CPU::op_ld_hl_sp_n(); d_dispatcher();
lf9: CPU::op_ld_sp_hl(); d_dispatcher();
lfa: CPU::op_ld_a_nn(); d_dispatcher();
lfb: CPU::op_ei(); d_dispatcher();
lfc: CPU::op_xx(); d_dispatcher();
lfd: CPU::op_xx(); d_dispatcher();
lfe: CPU::op_cp_a_n(); d_dispatcher();
lff: CPU::op_rst_n<0x38>(); d_dispatcher();
}

void CPU::Main() {
    cpu.exec();
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


// use test_for_interrupt macro instead 
inline void CPU::interrupt_test() {
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

inline void CPU::interrupt_exec(uint16 pc) {
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
  power_processor();

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
void CPU::power_processor() {
  r.halt = false;
  r.stop = false;
  r.ei = false;
  r.ime = false;
}
}



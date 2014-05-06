#include <processor/processor.hpp>
#include "lr35902.hpp"
#include <iostream>
#include <fstream>
#include <ostream>

namespace Processor {

#include "instructions.cpp"
#include "disassembler.cpp"
#include "serialization.cpp"

void LR35902::dump() {
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

void LR35902::power() {
  r.halt = false;
  r.stop = false;
  r.ei = false;
  r.ime = false;

  // normal instructions
  instructions[0x00] = &LR35902::op_nop;
  instructions[0x01] = &LR35902::op_ld_rr_nn<BC>;
  instructions[0x02] = &LR35902::op_ld_rr_a<BC>;
  instructions[0x03] = &LR35902::op_inc_rr<BC>;
  instructions[0x04] = &LR35902::op_inc_r<B>;
  instructions[0x05] = &LR35902::op_dec_r<B>;
  instructions[0x06] = &LR35902::op_ld_r_n<B>;
  instructions[0x07] = &LR35902::op_rlca;
  instructions[0x08] = &LR35902::op_ld_nn_sp;
  instructions[0x09] = &LR35902::op_add_hl_rr<BC>;
  instructions[0x0a] = &LR35902::op_ld_a_rr<BC>;
  instructions[0x0b] = &LR35902::op_dec_rr<BC>;
  instructions[0x0c] = &LR35902::op_inc_r<C>;
  instructions[0x0d] = &LR35902::op_dec_r<C>;
  instructions[0x0e] = &LR35902::op_ld_r_n<C>;
  instructions[0x0f] = &LR35902::op_rrca;
  instructions[0x10] = &LR35902::op_stop;
  instructions[0x11] = &LR35902::op_ld_rr_nn<DE>;
  instructions[0x12] = &LR35902::op_ld_rr_a<DE>;
  instructions[0x13] = &LR35902::op_inc_rr<DE>;
  instructions[0x14] = &LR35902::op_inc_r<D>;
  instructions[0x15] = &LR35902::op_dec_r<D>;
  instructions[0x16] = &LR35902::op_ld_r_n<D>;
  instructions[0x17] = &LR35902::op_rla;
  instructions[0x18] = &LR35902::op_jr_n;
  instructions[0x19] = &LR35902::op_add_hl_rr<DE>;
  instructions[0x1a] = &LR35902::op_ld_a_rr<DE>;
  instructions[0x1b] = &LR35902::op_dec_rr<DE>;
  instructions[0x1c] = &LR35902::op_inc_r<E>;
  instructions[0x1d] = &LR35902::op_dec_r<E>;
  instructions[0x1e] = &LR35902::op_ld_r_n<E>;
  instructions[0x1f] = &LR35902::op_rra;
  instructions[0x20] = &LR35902::op_jr_f_n<ZF, 0>;
  instructions[0x21] = &LR35902::op_ld_rr_nn<HL>;
  instructions[0x22] = &LR35902::op_ldi_hl_a;
  instructions[0x23] = &LR35902::op_inc_rr<HL>;
  instructions[0x24] = &LR35902::op_inc_r<H>;
  instructions[0x25] = &LR35902::op_dec_r<H>;
  instructions[0x26] = &LR35902::op_ld_r_n<H>;
  instructions[0x27] = &LR35902::op_daa;
  instructions[0x28] = &LR35902::op_jr_f_n<ZF, 1>;
  instructions[0x29] = &LR35902::op_add_hl_rr<HL>;
  instructions[0x2a] = &LR35902::op_ldi_a_hl;
  instructions[0x2b] = &LR35902::op_dec_rr<HL>;
  instructions[0x2c] = &LR35902::op_inc_r<L>;
  instructions[0x2d] = &LR35902::op_dec_r<L>;
  instructions[0x2e] = &LR35902::op_ld_r_n<L>;
  instructions[0x2f] = &LR35902::op_cpl;
  instructions[0x30] = &LR35902::op_jr_f_n<CF, 0>;
  instructions[0x31] = &LR35902::op_ld_rr_nn<SP>;
  instructions[0x32] = &LR35902::op_ldd_hl_a;
  instructions[0x33] = &LR35902::op_inc_rr<SP>;
  instructions[0x34] = &LR35902::op_inc_hl;
  instructions[0x35] = &LR35902::op_dec_hl;
  instructions[0x36] = &LR35902::op_ld_hl_n;
  instructions[0x37] = &LR35902::op_scf;
  instructions[0x38] = &LR35902::op_jr_f_n<CF, 1>;
  instructions[0x39] = &LR35902::op_add_hl_rr<SP>;
  instructions[0x3a] = &LR35902::op_ldd_a_hl;
  instructions[0x3b] = &LR35902::op_dec_rr<SP>;
  instructions[0x3c] = &LR35902::op_inc_r<A>;
  instructions[0x3d] = &LR35902::op_dec_r<A>;
  instructions[0x3e] = &LR35902::op_ld_r_n<A>;
  instructions[0x3f] = &LR35902::op_ccf;
  instructions[0x40] = &LR35902::op_ld_r_r<B, B>;
  instructions[0x41] = &LR35902::op_ld_r_r<B, C>;
  instructions[0x42] = &LR35902::op_ld_r_r<B, D>;
  instructions[0x43] = &LR35902::op_ld_r_r<B, E>;
  instructions[0x44] = &LR35902::op_ld_r_r<B, H>;
  instructions[0x45] = &LR35902::op_ld_r_r<B, L>;
  instructions[0x46] = &LR35902::op_ld_r_hl<B>;
  instructions[0x47] = &LR35902::op_ld_r_r<B, A>;
  instructions[0x48] = &LR35902::op_ld_r_r<C, B>;
  instructions[0x49] = &LR35902::op_ld_r_r<C, C>;
  instructions[0x4a] = &LR35902::op_ld_r_r<C, D>;
  instructions[0x4b] = &LR35902::op_ld_r_r<C, E>;
  instructions[0x4c] = &LR35902::op_ld_r_r<C, H>;
  instructions[0x4d] = &LR35902::op_ld_r_r<C, L>;
  instructions[0x4e] = &LR35902::op_ld_r_hl<C>;
  instructions[0x4f] = &LR35902::op_ld_r_r<C, A>;
  instructions[0x50] = &LR35902::op_ld_r_r<D, B>;
  instructions[0x51] = &LR35902::op_ld_r_r<D, C>;
  instructions[0x52] = &LR35902::op_ld_r_r<D, D>;
  instructions[0x53] = &LR35902::op_ld_r_r<D, E>;
  instructions[0x54] = &LR35902::op_ld_r_r<D, H>;
  instructions[0x55] = &LR35902::op_ld_r_r<D, L>;
  instructions[0x56] = &LR35902::op_ld_r_hl<D>;
  instructions[0x57] = &LR35902::op_ld_r_r<D, A>;
  instructions[0x58] = &LR35902::op_ld_r_r<E, B>;
  instructions[0x59] = &LR35902::op_ld_r_r<E, C>;
  instructions[0x5a] = &LR35902::op_ld_r_r<E, D>;
  instructions[0x5b] = &LR35902::op_ld_r_r<E, E>;
  instructions[0x5c] = &LR35902::op_ld_r_r<E, H>;
  instructions[0x5d] = &LR35902::op_ld_r_r<E, L>;
  instructions[0x5e] = &LR35902::op_ld_r_hl<E>;
  instructions[0x5f] = &LR35902::op_ld_r_r<E, A>;
  instructions[0x60] = &LR35902::op_ld_r_r<H, B>;
  instructions[0x61] = &LR35902::op_ld_r_r<H, C>;
  instructions[0x62] = &LR35902::op_ld_r_r<H, D>;
  instructions[0x63] = &LR35902::op_ld_r_r<H, E>;
  instructions[0x64] = &LR35902::op_ld_r_r<H, H>;
  instructions[0x65] = &LR35902::op_ld_r_r<H, L>;
  instructions[0x66] = &LR35902::op_ld_r_hl<H>;
  instructions[0x67] = &LR35902::op_ld_r_r<H, A>;
  instructions[0x68] = &LR35902::op_ld_r_r<L, B>;
  instructions[0x69] = &LR35902::op_ld_r_r<L, C>;
  instructions[0x6a] = &LR35902::op_ld_r_r<L, D>;
  instructions[0x6b] = &LR35902::op_ld_r_r<L, E>;
  instructions[0x6c] = &LR35902::op_ld_r_r<L, H>;
  instructions[0x6d] = &LR35902::op_ld_r_r<L, L>;
  instructions[0x6e] = &LR35902::op_ld_r_hl<L>;
  instructions[0x6f] = &LR35902::op_ld_r_r<L, A>;
  instructions[0x70] = &LR35902::op_ld_hl_r<B>;
  instructions[0x71] = &LR35902::op_ld_hl_r<C>;
  instructions[0x72] = &LR35902::op_ld_hl_r<D>;
  instructions[0x73] = &LR35902::op_ld_hl_r<E>;
  instructions[0x74] = &LR35902::op_ld_hl_r<H>;
  instructions[0x75] = &LR35902::op_ld_hl_r<L>;
  instructions[0x76] = &LR35902::op_halt;
  instructions[0x77] = &LR35902::op_ld_hl_r<A>;
  instructions[0x78] = &LR35902::op_ld_r_r<A, B>;
  instructions[0x79] = &LR35902::op_ld_r_r<A, C>;
  instructions[0x7a] = &LR35902::op_ld_r_r<A, D>;
  instructions[0x7b] = &LR35902::op_ld_r_r<A, E>;
  instructions[0x7c] = &LR35902::op_ld_r_r<A, H>;
  instructions[0x7d] = &LR35902::op_ld_r_r<A, L>;
  instructions[0x7e] = &LR35902::op_ld_r_hl<A>;
  instructions[0x7f] = &LR35902::op_ld_r_r<A, A>;
  instructions[0x80] = &LR35902::op_add_a_r<B>;
  instructions[0x81] = &LR35902::op_add_a_r<C>;
  instructions[0x82] = &LR35902::op_add_a_r<D>;
  instructions[0x83] = &LR35902::op_add_a_r<E>;
  instructions[0x84] = &LR35902::op_add_a_r<H>;
  instructions[0x85] = &LR35902::op_add_a_r<L>;
  instructions[0x86] = &LR35902::op_add_a_hl;
  instructions[0x87] = &LR35902::op_add_a_r<A>;
  instructions[0x88] = &LR35902::op_adc_a_r<B>;
  instructions[0x89] = &LR35902::op_adc_a_r<C>;
  instructions[0x8a] = &LR35902::op_adc_a_r<D>;
  instructions[0x8b] = &LR35902::op_adc_a_r<E>;
  instructions[0x8c] = &LR35902::op_adc_a_r<H>;
  instructions[0x8d] = &LR35902::op_adc_a_r<L>;
  instructions[0x8e] = &LR35902::op_adc_a_hl;
  instructions[0x8f] = &LR35902::op_adc_a_r<A>;
  instructions[0x90] = &LR35902::op_sub_a_r<B>;
  instructions[0x91] = &LR35902::op_sub_a_r<C>;
  instructions[0x92] = &LR35902::op_sub_a_r<D>;
  instructions[0x93] = &LR35902::op_sub_a_r<E>;
  instructions[0x94] = &LR35902::op_sub_a_r<H>;
  instructions[0x95] = &LR35902::op_sub_a_r<L>;
  instructions[0x96] = &LR35902::op_sub_a_hl;
  instructions[0x97] = &LR35902::op_sub_a_r<A>;
  instructions[0x98] = &LR35902::op_sbc_a_r<B>;
  instructions[0x99] = &LR35902::op_sbc_a_r<C>;
  instructions[0x9a] = &LR35902::op_sbc_a_r<D>;
  instructions[0x9b] = &LR35902::op_sbc_a_r<E>;
  instructions[0x9c] = &LR35902::op_sbc_a_r<H>;
  instructions[0x9d] = &LR35902::op_sbc_a_r<L>;
  instructions[0x9e] = &LR35902::op_sbc_a_hl;
  instructions[0x9f] = &LR35902::op_sbc_a_r<A>;
  instructions[0xa0] = &LR35902::op_and_a_r<B>;
  instructions[0xa1] = &LR35902::op_and_a_r<C>;
  instructions[0xa2] = &LR35902::op_and_a_r<D>;
  instructions[0xa3] = &LR35902::op_and_a_r<E>;
  instructions[0xa4] = &LR35902::op_and_a_r<H>;
  instructions[0xa5] = &LR35902::op_and_a_r<L>;
  instructions[0xa6] = &LR35902::op_and_a_hl;
  instructions[0xa7] = &LR35902::op_and_a_r<A>;
  instructions[0xa8] = &LR35902::op_xor_a_r<B>;
  instructions[0xa9] = &LR35902::op_xor_a_r<C>;
  instructions[0xaa] = &LR35902::op_xor_a_r<D>;
  instructions[0xab] = &LR35902::op_xor_a_r<E>;
  instructions[0xac] = &LR35902::op_xor_a_r<H>;
  instructions[0xad] = &LR35902::op_xor_a_r<L>;
  instructions[0xae] = &LR35902::op_xor_a_hl;
  instructions[0xaf] = &LR35902::op_xor_a_r<A>;
  instructions[0xb0] = &LR35902::op_or_a_r<B>;
  instructions[0xb1] = &LR35902::op_or_a_r<C>;
  instructions[0xb2] = &LR35902::op_or_a_r<D>;
  instructions[0xb3] = &LR35902::op_or_a_r<E>;
  instructions[0xb4] = &LR35902::op_or_a_r<H>;
  instructions[0xb5] = &LR35902::op_or_a_r<L>;
  instructions[0xb6] = &LR35902::op_or_a_hl;
  instructions[0xb7] = &LR35902::op_or_a_r<A>;
  instructions[0xb8] = &LR35902::op_cp_a_r<B>;
  instructions[0xb9] = &LR35902::op_cp_a_r<C>;
  instructions[0xba] = &LR35902::op_cp_a_r<D>;
  instructions[0xbb] = &LR35902::op_cp_a_r<E>;
  instructions[0xbc] = &LR35902::op_cp_a_r<H>;
  instructions[0xbd] = &LR35902::op_cp_a_r<L>;
  instructions[0xbe] = &LR35902::op_cp_a_hl;
  instructions[0xbf] = &LR35902::op_cp_a_r<A>;
  instructions[0xc0] = &LR35902::op_ret_f<ZF, 0>;
  instructions[0xc1] = &LR35902::op_pop_rr<BC>;
  instructions[0xc2] = &LR35902::op_jp_f_nn<ZF, 0>;
  instructions[0xc3] = &LR35902::op_jp_nn;
  instructions[0xc4] = &LR35902::op_call_f_nn<ZF, 0>;
  instructions[0xc5] = &LR35902::op_push_rr<BC>;
  instructions[0xc6] = &LR35902::op_add_a_n;
  instructions[0xc7] = &LR35902::op_rst_n<0x00>;
  instructions[0xc8] = &LR35902::op_ret_f<ZF, 1>;
  instructions[0xc9] = &LR35902::op_ret;
  instructions[0xca] = &LR35902::op_jp_f_nn<ZF, 1>;
  instructions[0xcb] = &LR35902::op_cb;
  instructions[0xcc] = &LR35902::op_call_f_nn<ZF, 1>;
  instructions[0xcd] = &LR35902::op_call_nn;
  instructions[0xce] = &LR35902::op_adc_a_n;
  instructions[0xcf] = &LR35902::op_rst_n<0x08>;
  instructions[0xd0] = &LR35902::op_ret_f<CF, 0>;
  instructions[0xd1] = &LR35902::op_pop_rr<DE>;
  instructions[0xd2] = &LR35902::op_jp_f_nn<CF, 0>;
  instructions[0xd3] = &LR35902::op_xx;
  instructions[0xd4] = &LR35902::op_call_f_nn<CF, 0>;
  instructions[0xd5] = &LR35902::op_push_rr<DE>;
  instructions[0xd6] = &LR35902::op_sub_a_n;
  instructions[0xd7] = &LR35902::op_rst_n<0x10>;
  instructions[0xd8] = &LR35902::op_ret_f<CF, 1>;
  instructions[0xd9] = &LR35902::op_reti;
  instructions[0xda] = &LR35902::op_jp_f_nn<CF, 1>;
  instructions[0xdb] = &LR35902::op_xx;
  instructions[0xdc] = &LR35902::op_call_f_nn<CF, 1>;
  instructions[0xdd] = &LR35902::op_xx;
  instructions[0xde] = &LR35902::op_sbc_a_n;
  instructions[0xdf] = &LR35902::op_rst_n<0x18>;
  instructions[0xe0] = &LR35902::op_ld_ffn_a;
  instructions[0xe1] = &LR35902::op_pop_rr<HL>;
  instructions[0xe2] = &LR35902::op_ld_ffc_a;
  instructions[0xe3] = &LR35902::op_xx;
  instructions[0xe4] = &LR35902::op_xx;
  instructions[0xe5] = &LR35902::op_push_rr<HL>;
  instructions[0xe6] = &LR35902::op_and_a_n;
  instructions[0xe7] = &LR35902::op_rst_n<0x20>;
  instructions[0xe8] = &LR35902::op_add_sp_n;
  instructions[0xe9] = &LR35902::op_jp_hl;
  instructions[0xea] = &LR35902::op_ld_nn_a;
  instructions[0xeb] = &LR35902::op_xx;
  instructions[0xec] = &LR35902::op_xx;
  instructions[0xed] = &LR35902::op_xx;
  instructions[0xee] = &LR35902::op_xor_a_n;
  instructions[0xef] = &LR35902::op_rst_n<0x28>;
  instructions[0xf0] = &LR35902::op_ld_a_ffn;
  instructions[0xf1] = &LR35902::op_pop_rr<AF>;
  instructions[0xf2] = &LR35902::op_ld_a_ffc;
  instructions[0xf3] = &LR35902::op_di;
  instructions[0xf4] = &LR35902::op_xx;
  instructions[0xf5] = &LR35902::op_push_rr<AF>;
  instructions[0xf6] = &LR35902::op_or_a_n;
  instructions[0xf7] = &LR35902::op_rst_n<0x30>;
  instructions[0xf8] = &LR35902::op_ld_hl_sp_n;
  instructions[0xf9] = &LR35902::op_ld_sp_hl;
  instructions[0xfa] = &LR35902::op_ld_a_nn;
  instructions[0xfb] = &LR35902::op_ei;
  instructions[0xfc] = &LR35902::op_xx;
  instructions[0xfd] = &LR35902::op_xx;
  instructions[0xfe] = &LR35902::op_cp_a_n;
  instructions[0xff] = &LR35902::op_rst_n<0x38>;

  // cb instructions
  cb_instructions[0x00] =&LR35902::op_rlc_r<B>;
  cb_instructions[0x01] =&LR35902::op_rlc_r<C>;
  cb_instructions[0x02] =&LR35902::op_rlc_r<D>;
  cb_instructions[0x03] =&LR35902::op_rlc_r<E>;
  cb_instructions[0x04] =&LR35902::op_rlc_r<H>;
  cb_instructions[0x05] =&LR35902::op_rlc_r<L>;
  cb_instructions[0x06] =&LR35902::op_rlc_hl;
  cb_instructions[0x07] =&LR35902::op_rlc_r<A>;
  cb_instructions[0x08] =&LR35902::op_rrc_r<B>;
  cb_instructions[0x09] =&LR35902::op_rrc_r<C>;
  cb_instructions[0x0a] =&LR35902::op_rrc_r<D>;
  cb_instructions[0x0b] =&LR35902::op_rrc_r<E>;
  cb_instructions[0x0c] =&LR35902::op_rrc_r<H>;
  cb_instructions[0x0d] =&LR35902::op_rrc_r<L>;
  cb_instructions[0x0e] =&LR35902::op_rrc_hl;
  cb_instructions[0x0f] =&LR35902::op_rrc_r<A>;
  cb_instructions[0x10] =&LR35902::op_rl_r<B>;
  cb_instructions[0x11] =&LR35902::op_rl_r<C>;
  cb_instructions[0x12] =&LR35902::op_rl_r<D>;
  cb_instructions[0x13] =&LR35902::op_rl_r<E>;
  cb_instructions[0x14] =&LR35902::op_rl_r<H>;
  cb_instructions[0x15] =&LR35902::op_rl_r<L>;
  cb_instructions[0x16] =&LR35902::op_rl_hl;
  cb_instructions[0x17] =&LR35902::op_rl_r<A>;
  cb_instructions[0x18] =&LR35902::op_rr_r<B>;
  cb_instructions[0x19] =&LR35902::op_rr_r<C>;
  cb_instructions[0x1a] =&LR35902::op_rr_r<D>;
  cb_instructions[0x1b] =&LR35902::op_rr_r<E>;
  cb_instructions[0x1c] =&LR35902::op_rr_r<H>;
  cb_instructions[0x1d] =&LR35902::op_rr_r<L>;
  cb_instructions[0x1e] =&LR35902::op_rr_hl;
  cb_instructions[0x1f] =&LR35902::op_rr_r<A>;
  cb_instructions[0x20] =&LR35902::op_sla_r<B>;
  cb_instructions[0x21] =&LR35902::op_sla_r<C>;
  cb_instructions[0x22] =&LR35902::op_sla_r<D>;
  cb_instructions[0x23] =&LR35902::op_sla_r<E>;
  cb_instructions[0x24] =&LR35902::op_sla_r<H>;
  cb_instructions[0x25] =&LR35902::op_sla_r<L>;
  cb_instructions[0x26] =&LR35902::op_sla_hl;
  cb_instructions[0x27] =&LR35902::op_sla_r<A>;
  cb_instructions[0x28] =&LR35902::op_sra_r<B>;
  cb_instructions[0x29] =&LR35902::op_sra_r<C>;
  cb_instructions[0x2a] =&LR35902::op_sra_r<D>;
  cb_instructions[0x2b] =&LR35902::op_sra_r<E>;
  cb_instructions[0x2c] =&LR35902::op_sra_r<H>;
  cb_instructions[0x2d] =&LR35902::op_sra_r<L>;
  cb_instructions[0x2e] =&LR35902::op_sra_hl;
  cb_instructions[0x2f] =&LR35902::op_sra_r<A>;
  cb_instructions[0x30] =&LR35902::op_swap_r<B>;
  cb_instructions[0x31] =&LR35902::op_swap_r<C>;
  cb_instructions[0x32] =&LR35902::op_swap_r<D>;
  cb_instructions[0x33] =&LR35902::op_swap_r<E>;
  cb_instructions[0x34] =&LR35902::op_swap_r<H>;
  cb_instructions[0x35] =&LR35902::op_swap_r<L>;
  cb_instructions[0x36] =&LR35902::op_swap_hl;
  cb_instructions[0x37] =&LR35902::op_swap_r<A>;
  cb_instructions[0x38] =&LR35902::op_srl_r<B>;
  cb_instructions[0x39] =&LR35902::op_srl_r<C>;
  cb_instructions[0x3a] =&LR35902::op_srl_r<D>;
  cb_instructions[0x3b] =&LR35902::op_srl_r<E>;
  cb_instructions[0x3c] =&LR35902::op_srl_r<H>;
  cb_instructions[0x3d] =&LR35902::op_srl_r<L>;
  cb_instructions[0x3e] =&LR35902::op_srl_hl;
  cb_instructions[0x3f] =&LR35902::op_srl_r<A>;
  cb_instructions[0x40] =&LR35902::op_bit_n_r<0, B>;
  cb_instructions[0x41] =&LR35902::op_bit_n_r<0, C>;
  cb_instructions[0x42] =&LR35902::op_bit_n_r<0, D>;
  cb_instructions[0x43] =&LR35902::op_bit_n_r<0, E>;
  cb_instructions[0x44] =&LR35902::op_bit_n_r<0, H>;
  cb_instructions[0x45] =&LR35902::op_bit_n_r<0, L>;
  cb_instructions[0x46] =&LR35902::op_bit_n_hl<0>;
  cb_instructions[0x47] =&LR35902::op_bit_n_r<0, A>;
  cb_instructions[0x48] =&LR35902::op_bit_n_r<1, B>;
  cb_instructions[0x49] =&LR35902::op_bit_n_r<1, C>;
  cb_instructions[0x4a] =&LR35902::op_bit_n_r<1, D>;
  cb_instructions[0x4b] =&LR35902::op_bit_n_r<1, E>;
  cb_instructions[0x4c] =&LR35902::op_bit_n_r<1, H>;
  cb_instructions[0x4d] =&LR35902::op_bit_n_r<1, L>;
  cb_instructions[0x4e] =&LR35902::op_bit_n_hl<1>;
  cb_instructions[0x4f] =&LR35902::op_bit_n_r<1, A>;
  cb_instructions[0x50] =&LR35902::op_bit_n_r<2, B>;
  cb_instructions[0x51] =&LR35902::op_bit_n_r<2, C>;
  cb_instructions[0x52] =&LR35902::op_bit_n_r<2, D>;
  cb_instructions[0x53] =&LR35902::op_bit_n_r<2, E>;
  cb_instructions[0x54] =&LR35902::op_bit_n_r<2, H>;
  cb_instructions[0x55] =&LR35902::op_bit_n_r<2, L>;
  cb_instructions[0x56] =&LR35902::op_bit_n_hl<2>;
  cb_instructions[0x57] =&LR35902::op_bit_n_r<2, A>;
  cb_instructions[0x58] =&LR35902::op_bit_n_r<3, B>;
  cb_instructions[0x59] =&LR35902::op_bit_n_r<3, C>;
  cb_instructions[0x5a] =&LR35902::op_bit_n_r<3, D>;
  cb_instructions[0x5b] =&LR35902::op_bit_n_r<3, E>;
  cb_instructions[0x5c] =&LR35902::op_bit_n_r<3, H>;
  cb_instructions[0x5d] =&LR35902::op_bit_n_r<3, L>;
  cb_instructions[0x5e] =&LR35902::op_bit_n_hl<3>;
  cb_instructions[0x5f] =&LR35902::op_bit_n_r<3, A>;
  cb_instructions[0x60] =&LR35902::op_bit_n_r<4, B>;
  cb_instructions[0x61] =&LR35902::op_bit_n_r<4, C>;
  cb_instructions[0x62] =&LR35902::op_bit_n_r<4, D>;
  cb_instructions[0x63] =&LR35902::op_bit_n_r<4, E>;
  cb_instructions[0x64] =&LR35902::op_bit_n_r<4, H>;
  cb_instructions[0x65] =&LR35902::op_bit_n_r<4, L>;
  cb_instructions[0x66] =&LR35902::op_bit_n_hl<4>;
  cb_instructions[0x67] =&LR35902::op_bit_n_r<4, A>;
  cb_instructions[0x68] =&LR35902::op_bit_n_r<5, B>;
  cb_instructions[0x69] =&LR35902::op_bit_n_r<5, C>;
  cb_instructions[0x6a] =&LR35902::op_bit_n_r<5, D>;
  cb_instructions[0x6b] =&LR35902::op_bit_n_r<5, E>;
  cb_instructions[0x6c] =&LR35902::op_bit_n_r<5, H>;
  cb_instructions[0x6d] =&LR35902::op_bit_n_r<5, L>;
  cb_instructions[0x6e] =&LR35902::op_bit_n_hl<5>;
  cb_instructions[0x6f] =&LR35902::op_bit_n_r<5, A>;
  cb_instructions[0x70] =&LR35902::op_bit_n_r<6, B>;
  cb_instructions[0x71] =&LR35902::op_bit_n_r<6, C>;
  cb_instructions[0x72] =&LR35902::op_bit_n_r<6, D>;
  cb_instructions[0x73] =&LR35902::op_bit_n_r<6, E>;
  cb_instructions[0x74] =&LR35902::op_bit_n_r<6, H>;
  cb_instructions[0x75] =&LR35902::op_bit_n_r<6, L>;
  cb_instructions[0x76] =&LR35902::op_bit_n_hl<6>;
  cb_instructions[0x77] =&LR35902::op_bit_n_r<6, A>;
  cb_instructions[0x78] =&LR35902::op_bit_n_r<7, B>;
  cb_instructions[0x79] =&LR35902::op_bit_n_r<7, C>;
  cb_instructions[0x7a] =&LR35902::op_bit_n_r<7, D>;
  cb_instructions[0x7b] =&LR35902::op_bit_n_r<7, E>;
  cb_instructions[0x7c] =&LR35902::op_bit_n_r<7, H>;
  cb_instructions[0x7d] =&LR35902::op_bit_n_r<7, L>;
  cb_instructions[0x7e] =&LR35902::op_bit_n_hl<7>;
  cb_instructions[0x7f] =&LR35902::op_bit_n_r<7, A>;
  cb_instructions[0x80] =&LR35902::op_res_n_r<0, B>;
  cb_instructions[0x81] =&LR35902::op_res_n_r<0, C>;
  cb_instructions[0x82] =&LR35902::op_res_n_r<0, D>;
  cb_instructions[0x83] =&LR35902::op_res_n_r<0, E>;
  cb_instructions[0x84] =&LR35902::op_res_n_r<0, H>;
  cb_instructions[0x85] =&LR35902::op_res_n_r<0, L>;
  cb_instructions[0x86] =&LR35902::op_res_n_hl<0>;
  cb_instructions[0x87] =&LR35902::op_res_n_r<0, A>;
  cb_instructions[0x88] =&LR35902::op_res_n_r<1, B>;
  cb_instructions[0x89] =&LR35902::op_res_n_r<1, C>;
  cb_instructions[0x8a] =&LR35902::op_res_n_r<1, D>;
  cb_instructions[0x8b] =&LR35902::op_res_n_r<1, E>;
  cb_instructions[0x8c] =&LR35902::op_res_n_r<1, H>;
  cb_instructions[0x8d] =&LR35902::op_res_n_r<1, L>;
  cb_instructions[0x8e] =&LR35902::op_res_n_hl<1>;
  cb_instructions[0x8f] =&LR35902::op_res_n_r<1, A>;
  cb_instructions[0x90] =&LR35902::op_res_n_r<2, B>;
  cb_instructions[0x91] =&LR35902::op_res_n_r<2, C>;
  cb_instructions[0x92] =&LR35902::op_res_n_r<2, D>;
  cb_instructions[0x93] =&LR35902::op_res_n_r<2, E>;
  cb_instructions[0x94] =&LR35902::op_res_n_r<2, H>;
  cb_instructions[0x95] =&LR35902::op_res_n_r<2, L>;
  cb_instructions[0x96] =&LR35902::op_res_n_hl<2>;
  cb_instructions[0x97] =&LR35902::op_res_n_r<2, A>;
  cb_instructions[0x98] =&LR35902::op_res_n_r<3, B>;
  cb_instructions[0x99] =&LR35902::op_res_n_r<3, C>;
  cb_instructions[0x9a] =&LR35902::op_res_n_r<3, D>;
  cb_instructions[0x9b] =&LR35902::op_res_n_r<3, E>;
  cb_instructions[0x9c] =&LR35902::op_res_n_r<3, H>;
  cb_instructions[0x9d] =&LR35902::op_res_n_r<3, L>;
  cb_instructions[0x9e] =&LR35902::op_res_n_hl<3>;
  cb_instructions[0x9f] =&LR35902::op_res_n_r<3, A>;
  cb_instructions[0xa0] =&LR35902::op_res_n_r<4, B>;
  cb_instructions[0xa1] =&LR35902::op_res_n_r<4, C>;
  cb_instructions[0xa2] =&LR35902::op_res_n_r<4, D>;
  cb_instructions[0xa3] =&LR35902::op_res_n_r<4, E>;
  cb_instructions[0xa4] =&LR35902::op_res_n_r<4, H>;
  cb_instructions[0xa5] =&LR35902::op_res_n_r<4, L>;
  cb_instructions[0xa6] =&LR35902::op_res_n_hl<4>;
  cb_instructions[0xa7] =&LR35902::op_res_n_r<4, A>;
  cb_instructions[0xa8] =&LR35902::op_res_n_r<5, B>;
  cb_instructions[0xa9] =&LR35902::op_res_n_r<5, C>;
  cb_instructions[0xaa] =&LR35902::op_res_n_r<5, D>;
  cb_instructions[0xab] =&LR35902::op_res_n_r<5, E>;
  cb_instructions[0xac] =&LR35902::op_res_n_r<5, H>;
  cb_instructions[0xad] =&LR35902::op_res_n_r<5, L>;
  cb_instructions[0xae] =&LR35902::op_res_n_hl<5>;
  cb_instructions[0xaf] =&LR35902::op_res_n_r<5, A>;
  cb_instructions[0xb0] =&LR35902::op_res_n_r<6, B>;
  cb_instructions[0xb1] =&LR35902::op_res_n_r<6, C>;
  cb_instructions[0xb2] =&LR35902::op_res_n_r<6, D>;
  cb_instructions[0xb3] =&LR35902::op_res_n_r<6, E>;
  cb_instructions[0xb4] =&LR35902::op_res_n_r<6, H>;
  cb_instructions[0xb5] =&LR35902::op_res_n_r<6, L>;
  cb_instructions[0xb6] =&LR35902::op_res_n_hl<6>;
  cb_instructions[0xb7] =&LR35902::op_res_n_r<6, A>;
  cb_instructions[0xb8] =&LR35902::op_res_n_r<7, B>;
  cb_instructions[0xb9] =&LR35902::op_res_n_r<7, C>;
  cb_instructions[0xba] =&LR35902::op_res_n_r<7, D>;
  cb_instructions[0xbb] =&LR35902::op_res_n_r<7, E>;
  cb_instructions[0xbc] =&LR35902::op_res_n_r<7, H>;
  cb_instructions[0xbd] =&LR35902::op_res_n_r<7, L>;
  cb_instructions[0xbe] =&LR35902::op_res_n_hl<7>;
  cb_instructions[0xbf] =&LR35902::op_res_n_r<7, A>;
  cb_instructions[0xc0] =&LR35902::op_set_n_r<0, B>;
  cb_instructions[0xc1] =&LR35902::op_set_n_r<0, C>;
  cb_instructions[0xc2] =&LR35902::op_set_n_r<0, D>;
  cb_instructions[0xc3] =&LR35902::op_set_n_r<0, E>;
  cb_instructions[0xc4] =&LR35902::op_set_n_r<0, H>;
  cb_instructions[0xc5] =&LR35902::op_set_n_r<0, L>;
  cb_instructions[0xc6] =&LR35902::op_set_n_hl<0>;
  cb_instructions[0xc7] =&LR35902::op_set_n_r<0, A>;
  cb_instructions[0xc8] =&LR35902::op_set_n_r<1, B>;
  cb_instructions[0xc9] =&LR35902::op_set_n_r<1, C>;
  cb_instructions[0xca] =&LR35902::op_set_n_r<1, D>;
  cb_instructions[0xcb] =&LR35902::op_set_n_r<1, E>;
  cb_instructions[0xcc] =&LR35902::op_set_n_r<1, H>;
  cb_instructions[0xcd] =&LR35902::op_set_n_r<1, L>;
  cb_instructions[0xce] =&LR35902::op_set_n_hl<1>;
  cb_instructions[0xcf] =&LR35902::op_set_n_r<1, A>;
  cb_instructions[0xd0] =&LR35902::op_set_n_r<2, B>;
  cb_instructions[0xd1] =&LR35902::op_set_n_r<2, C>;
  cb_instructions[0xd2] =&LR35902::op_set_n_r<2, D>;
  cb_instructions[0xd3] =&LR35902::op_set_n_r<2, E>;
  cb_instructions[0xd4] =&LR35902::op_set_n_r<2, H>;
  cb_instructions[0xd5] =&LR35902::op_set_n_r<2, L>;
  cb_instructions[0xd6] =&LR35902::op_set_n_hl<2>;
  cb_instructions[0xd7] =&LR35902::op_set_n_r<2, A>;
  cb_instructions[0xd8] =&LR35902::op_set_n_r<3, B>;
  cb_instructions[0xd9] =&LR35902::op_set_n_r<3, C>;
  cb_instructions[0xda] =&LR35902::op_set_n_r<3, D>;
  cb_instructions[0xdb] =&LR35902::op_set_n_r<3, E>;
  cb_instructions[0xdc] =&LR35902::op_set_n_r<3, H>;
  cb_instructions[0xdd] =&LR35902::op_set_n_r<3, L>;
  cb_instructions[0xde] =&LR35902::op_set_n_hl<3>;
  cb_instructions[0xdf] =&LR35902::op_set_n_r<3, A>;
  cb_instructions[0xe0] =&LR35902::op_set_n_r<4, B>;
  cb_instructions[0xe1] =&LR35902::op_set_n_r<4, C>;
  cb_instructions[0xe2] =&LR35902::op_set_n_r<4, D>;
  cb_instructions[0xe3] =&LR35902::op_set_n_r<4, E>;
  cb_instructions[0xe4] =&LR35902::op_set_n_r<4, H>;
  cb_instructions[0xe5] =&LR35902::op_set_n_r<4, L>;
  cb_instructions[0xe6] =&LR35902::op_set_n_hl<4>;
  cb_instructions[0xe7] =&LR35902::op_set_n_r<4, A>;
  cb_instructions[0xe8] =&LR35902::op_set_n_r<5, B>;
  cb_instructions[0xe9] =&LR35902::op_set_n_r<5, C>;
  cb_instructions[0xea] =&LR35902::op_set_n_r<5, D>;
  cb_instructions[0xeb] =&LR35902::op_set_n_r<5, E>;
  cb_instructions[0xec] =&LR35902::op_set_n_r<5, H>;
  cb_instructions[0xed] =&LR35902::op_set_n_r<5, L>;
  cb_instructions[0xee] =&LR35902::op_set_n_hl<5>;
  cb_instructions[0xef] =&LR35902::op_set_n_r<5, A>;
  cb_instructions[0xf0] =&LR35902::op_set_n_r<6, B>;
  cb_instructions[0xf1] =&LR35902::op_set_n_r<6, C>;
  cb_instructions[0xf2] =&LR35902::op_set_n_r<6, D>;
  cb_instructions[0xf3] =&LR35902::op_set_n_r<6, E>;
  cb_instructions[0xf4] =&LR35902::op_set_n_r<6, H>;
  cb_instructions[0xf5] =&LR35902::op_set_n_r<6, L>;
  cb_instructions[0xf6] =&LR35902::op_set_n_hl<6>;
  cb_instructions[0xf7] =&LR35902::op_set_n_r<6, A>;
  cb_instructions[0xf8] =&LR35902::op_set_n_r<7, B>;
  cb_instructions[0xf9] =&LR35902::op_set_n_r<7, C>;
  cb_instructions[0xfa] =&LR35902::op_set_n_r<7, D>;
  cb_instructions[0xfb] =&LR35902::op_set_n_r<7, E>;
  cb_instructions[0xfc] =&LR35902::op_set_n_r<7, H>;
  cb_instructions[0xfd] =&LR35902::op_set_n_r<7, L>;
  cb_instructions[0xfe] =&LR35902::op_set_n_hl<7>;
  cb_instructions[0xff] =&LR35902::op_set_n_r<7, A>;
}

void LR35902::exec() {
  uint8 opcode = op_read(r[PC]++);
  last_inst = opcode;
  inst_counter[opcode]++;
  instruction_count++;
  (this->*instructions[opcode])();
}

void LR35902::exec_cb() {
  cb_operation = true;
  uint8 opcode = op_read(r[PC]++);
  last_inst = opcode;
  cb_inst_counter[opcode]++;
  (this->*cb_instructions[opcode])();
}

}

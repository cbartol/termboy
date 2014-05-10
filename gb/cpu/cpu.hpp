struct CPU : Thread, MMIO {
    enum {
        A, F, AF,
        B, C, BC,
        D, E, DE,
        H, L, HL,
        SP, PC,
    };

    enum {
        ZF, NF, HF, CF,
    };

    //register base class
    //the idea here is to have all registers derive from a single base class.
    //this allows construction of opcodes that can take any register as input or output,
    //despite the fact that behind-the-scenes, special handling is done for eg: F, AF, HL, etc.
    //registers can also be chained together: eg af = 0x0000 writes both a and f.
    struct Register {
        virtual operator unsigned() const = 0;
        virtual unsigned operator=(unsigned x) = 0;
        Register& operator=(const Register& x) { operator=((unsigned)x); return *this; }

        unsigned operator++(int) { unsigned r = *this; operator=(*this + 1); return r; }
        unsigned operator--(int) { unsigned r = *this; operator=(*this - 1); return r; }
        unsigned operator++() { return operator=(*this + 1); }
        unsigned operator--() { return operator=(*this - 1); }

        unsigned operator |=(unsigned x) { return operator=(*this | x); }
        unsigned operator ^=(unsigned x) { return operator=(*this ^ x); }
        unsigned operator &=(unsigned x) { return operator=(*this & x); }

        unsigned operator<<=(unsigned x) { return operator=(*this << x); }
        unsigned operator>>=(unsigned x) { return operator=(*this >> x); }

        unsigned operator +=(unsigned x) { return operator=(*this + x); }
        unsigned operator -=(unsigned x) { return operator=(*this - x); }
        unsigned operator *=(unsigned x) { return operator=(*this * x); }
        unsigned operator /=(unsigned x) { return operator=(*this / x); }
        unsigned operator %=(unsigned x) { return operator=(*this % x); }
    };

    struct Register8 : Register {
        uint8 data;
        operator unsigned() const { return data; }
        unsigned operator=(unsigned x) { return data = x; }
    };

    struct RegisterF : Register {
        bool z, n, h, c;
        operator unsigned() const { return (z << 7) | (n << 6) | (h << 5) | (c << 4); }
        unsigned operator=(unsigned x) { z = x & 0x80; n = x & 0x40; h = x & 0x20; c = x & 0x10; return *this; }
        bool& operator[](unsigned r) {
            static bool* table[] = {&z, &n, &h, &c};
            return *table[r];
        }
    };

    struct Register16 : Register {
        uint16 data;
        operator unsigned() const { return data; }
        unsigned operator=(unsigned x) { return data = x; }
    };

    struct RegisterAF : Register {
        Register8& hi;
        RegisterF& lo;
        operator unsigned() const { return (hi << 8) | (lo << 0); }
        unsigned operator=(unsigned x) { hi = x >> 8; lo = x >> 0; return *this; }
        RegisterAF(Register8& hi, RegisterF& lo) : hi(hi), lo(lo) {}
    };

    struct RegisterW : Register {
        Register8& hi;
        Register8& lo;
        operator unsigned() const { return (hi << 8) | (lo << 0); }
        unsigned operator=(unsigned x) { hi = x >> 8; lo = x >> 0; return *this; }
        RegisterW(Register8& hi, Register8& lo) : hi(hi), lo(lo) {}
    };

    struct Registers {
        Register8  a;
        RegisterF  f;
        RegisterAF af;
        Register8  b;
        Register8  c;
        RegisterW  bc;
        Register8  d;
        Register8  e;
        RegisterW  de;
        Register8  h;
        Register8  l;
        RegisterW  hl;
        Register16 sp;
        Register16 pc;

        bool halt;
        bool stop;
        bool ei;
        bool ime;

        Register& operator[](unsigned r) {
            static Register* table[] = {&a, &f, &af, &b, &c, &bc, &d, &e, &de, &h, &l, &hl, &sp, &pc};
            return *table[r];
        }

        Registers() : af(a, f), bc(b, c), de(d, e), hl(h, l) {}
    } r;


    enum class Interrupt : unsigned {
        Vblank,
        Stat,
        Timer,
        Serial,
        Joypad,
    };

    void (CPU::*instructions[256])();
    void (CPU::*cb_instructions[256])(); 

    double global_time;
    double time = 0;
    double max_time;
    double cb_time;
    double cb_max_time;
    int instruction_count = 0;
    bool cb_operation = false;
    double synch_time = 0;


    void dump();
    void power_processor();
    void exec();
    void exec_cb();

    void get_times(double final_time, double synch);


    struct Status {
        unsigned clock;

        //$ff00  JOYP
        bool p15;
        bool p14;
        uint8 joyp;
        uint8 mlt_req;

        //$ff01  SB
        uint8 serial_data;
        unsigned serial_bits;

        //$ff02  SC
        bool serial_transfer;
        bool serial_clock;

        //$ff04  DIV
        uint8 div;

        //$ff05  TIMA
        uint8 tima;

        //$ff06  TMA
        uint8 tma;

        //$ff07  TAC
        bool timer_enable;
        unsigned timer_clock;

        //$ff0f  IF
        bool interrupt_request_joypad;
        bool interrupt_request_serial;
        bool interrupt_request_timer;
        bool interrupt_request_stat;
        bool interrupt_request_vblank;

        //$ff4d  KEY1
        bool speed_double;
        bool speed_switch;

        //$ff51,$ff52  HDMA1,HDMA2
        uint16 dma_source;

        //$ff53,$ff54  HDMA3,HDMA4
        uint16 dma_target;

        //$ff55  HDMA5
        bool dma_mode;
        uint16 dma_length;

        //$ff6c  ???
        uint8 ff6c;

        //$ff70  SVBK
        uint3 wram_bank;

        //$ff72-$ff75  ???
        uint8 ff72;
        uint8 ff73;
        uint8 ff74;
        uint8 ff75;

        //$ffff  IE
        bool interrupt_enable_joypad;
        bool interrupt_enable_serial;
        bool interrupt_enable_timer;
        bool interrupt_enable_stat;
        bool interrupt_enable_vblank;
    } status;

    uint8 wram[32768];  //GB=8192, GBC=32768
    uint8 hram[128];

    static void Main();
    void interrupt_raise(Interrupt id);
    void interrupt_test();
    void interrupt_exec(uint16 pc);
    bool stop();
    void power();
    void serialize(serializer&);

    //mmio.cpp
    unsigned wram_addr(uint16 addr) const;
    void mmio_joyp_poll();
    uint8 mmio_read(uint16 addr);
    void mmio_write(uint16 addr, uint8 data);

    //memory.cpp
    void op_io();
    uint8 op_read(uint16 addr);
    void op_write(uint16 addr, uint8 data);
    void cycle_edge();
    uint8 debugger_read(uint16 addr);

    //timing.cpp
    void add_clocks(unsigned clocks);
    void timer_262144hz();
    void timer_65536hz();
    void timer_16384hz();
    void timer_8192hz();
    void timer_4096hz();
    void hblank();

privileged:
    void op_xx();
    void op_cb();

    //8-bit load commands
    template<unsigned x, unsigned y> void op_ld_r_r();
    template<unsigned x> void op_ld_r_n();
    template<unsigned x> void op_ld_r_hl();
    template<unsigned x> void op_ld_hl_r();
    void op_ld_hl_n();
    template<unsigned x> void op_ld_a_rr();
    void op_ld_a_nn();
    template<unsigned x> void op_ld_rr_a();
    void op_ld_nn_a();
    void op_ld_a_ffn();
    void op_ld_ffn_a();
    void op_ld_a_ffc();
    void op_ld_ffc_a();
    void op_ldi_hl_a();
    void op_ldi_a_hl();
    void op_ldd_hl_a();
    void op_ldd_a_hl();

    //16-bit load commands
    template<unsigned x> void op_ld_rr_nn();
    void op_ld_nn_sp();
    void op_ld_sp_hl();
    template<unsigned x> void op_push_rr();
    template<unsigned x> void op_pop_rr();

    //8-bit arithmetic commands
    void opi_add_a(uint8 x);
    template<unsigned x> void op_add_a_r();
    void op_add_a_n();
    void op_add_a_hl();

    void opi_adc_a(uint8 x);
    template<unsigned x> void op_adc_a_r();
    void op_adc_a_n();
    void op_adc_a_hl();

    void opi_sub_a(uint8 x);
    template<unsigned x> void op_sub_a_r();
    void op_sub_a_n();
    void op_sub_a_hl();

    void opi_sbc_a(uint8 x);
    template<unsigned x> void op_sbc_a_r();
    void op_sbc_a_n();
    void op_sbc_a_hl();

    void opi_and_a(uint8 x);
    template<unsigned x> void op_and_a_r();
    void op_and_a_n();
    void op_and_a_hl();

    void opi_xor_a(uint8 x);
    template<unsigned x> void op_xor_a_r();
    void op_xor_a_n();
    void op_xor_a_hl();

    void opi_or_a(uint8 x);
    template<unsigned x> void op_or_a_r();
    void op_or_a_n();
    void op_or_a_hl();

    void opi_cp_a(uint8 x);
    template<unsigned x> void op_cp_a_r();
    void op_cp_a_n();
    void op_cp_a_hl();

    template<unsigned x> void op_inc_r();
    void op_inc_hl();
    template<unsigned x> void op_dec_r();
    void op_dec_hl();
    void op_daa();
    void op_cpl();

    //16-bit arithmetic commands
    template<unsigned x> void op_add_hl_rr();
    template<unsigned x> void op_inc_rr();
    template<unsigned x> void op_dec_rr();
    void op_add_sp_n();
    void op_ld_hl_sp_n();

    //rotate/shift commands
    void op_rlca();
    void op_rla();
    void op_rrca();
    void op_rra();
    template<unsigned x> void op_rlc_r();
    void op_rlc_hl();
    template<unsigned x> void op_rl_r();
    void op_rl_hl();
    template<unsigned x> void op_rrc_r();
    void op_rrc_hl();
    template<unsigned x> void op_rr_r();
    void op_rr_hl();
    template<unsigned x> void op_sla_r();
    void op_sla_hl();
    template<unsigned x> void op_swap_r();
    void op_swap_hl();
    template<unsigned x> void op_sra_r();
    void op_sra_hl();
    template<unsigned x> void op_srl_r();
    void op_srl_hl();

    //single-bit commands
    template<unsigned b, unsigned x> void op_bit_n_r();
    template<unsigned b> void op_bit_n_hl();
    template<unsigned b, unsigned x> void op_set_n_r();
    template<unsigned b> void op_set_n_hl();
    template<unsigned b, unsigned x> void op_res_n_r();
    template<unsigned b> void op_res_n_hl();

    //control commands
    void op_ccf();
    void op_scf();
    void op_nop();
    void op_halt();
    void op_stop();
    void op_di();
    void op_ei();

    //jump commands
    void op_jp_nn();
    void op_jp_hl();
    template<unsigned x, bool y> void op_jp_f_nn();
    void op_jr_n();
    template<unsigned x, bool y> void op_jr_f_n();
    void op_call_nn();
    template<unsigned x, bool y> void op_call_f_nn();
    void op_ret();
    template<unsigned x, bool y> void op_ret_f();
    void op_reti();
    template<unsigned n> void op_rst_n();
};

extern CPU cpu;

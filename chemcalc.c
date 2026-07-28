/*
 * chemcalc.c - a small CLI tool for chemistry calculations.
 *
 * Features:
 *   mw     - compute molecular weight (and mmol from a given mass)
 *   react  - interactively enter reactants + product(s), compute theoretical yield
 *   list   - list saved reactions
 *   show   - show details of one saved reaction
 *   yield  - record an actual yield and compute % yield
 *   rename - give a saved reaction a new (named) id
 *   delete - remove a saved reaction
 *
 * A reaction id may be a name you choose ("aldol-1") instead of the random
 * hex id; `react` asks for one and `rename` changes it afterwards.
 *
 * Formula input is case-insensitive: "C6H6", "c6h6" and "ccl4" all work.
 * Reactants and products may carry a leading coefficient (equivalents),
 * e.g. "2NaHCO3".  Reactions may have more than one product.
 *
 * Masses are handled in milligrams and amounts in millimoles throughout;
 * other units (g/kg/ug) are accepted on input and converted immediately.
 *
 * Data is stored as a simple tab-delimited text file at:
 *   ~/.local/share/chemcalc/reactions.dat
 *
 * Build:
 *   gcc -O2 -Wall -o chemcalc chemcalc.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

/* -------------------------------------------------------------------- */
/* Periodic table (atomic weights, g/mol)                               */
/* -------------------------------------------------------------------- */
typedef struct { const char *symbol; double weight; } Element;

static const Element PERIODIC_TABLE[] = {
    {"H",1.008},{"He",4.0026},{"Li",6.94},{"Be",9.0122},{"B",10.81},
    {"C",12.011},{"N",14.007},{"O",15.999},{"F",18.998},{"Ne",20.180},
    {"Na",22.990},{"Mg",24.305},{"Al",26.982},{"Si",28.085},{"P",30.974},
    {"S",32.06},{"Cl",35.45},{"Ar",39.948},{"K",39.098},{"Ca",40.078},
    {"Sc",44.956},{"Ti",47.867},{"V",50.942},{"Cr",51.996},{"Mn",54.938},
    {"Fe",55.845},{"Co",58.933},{"Ni",58.693},{"Cu",63.546},{"Zn",65.38},
    {"Ga",69.723},{"Ge",72.630},{"As",74.922},{"Se",78.971},{"Br",79.904},
    {"Kr",83.798},{"Rb",85.468},{"Sr",87.62},{"Y",88.906},{"Zr",91.224},
    {"Nb",92.906},{"Mo",95.95},{"Tc",98.0},{"Ru",101.07},{"Rh",102.91},
    {"Pd",106.42},{"Ag",107.87},{"Cd",112.41},{"In",114.82},{"Sn",118.71},
    {"Sb",121.76},{"Te",127.60},{"I",126.90},{"Xe",131.29},{"Cs",132.91},
    {"Ba",137.33},{"La",138.91},{"Ce",140.12},{"Pr",140.91},{"Nd",144.24},
    {"Pm",145.0},{"Sm",150.36},{"Eu",151.96},{"Gd",157.25},{"Tb",158.93},
    {"Dy",162.50},{"Ho",164.93},{"Er",167.26},{"Tm",168.93},{"Yb",173.05},
    {"Lu",174.97},{"Hf",178.49},{"Ta",180.95},{"W",183.84},{"Re",186.21},
    {"Os",190.23},{"Ir",192.22},{"Pt",195.08},{"Au",196.97},{"Hg",200.59},
    {"Tl",204.38},{"Pb",207.2},{"Bi",208.98},{"Po",209.0},{"At",210.0},
    {"Rn",222.0},{"Fr",223.0},{"Ra",226.0},{"Ac",227.0},{"Th",232.04},
    {"Pa",231.04},{"U",238.03},{"Np",237.0},{"Pu",244.0},
};
static const int PERIODIC_TABLE_LEN = sizeof(PERIODIC_TABLE) / sizeof(Element);

static double lookup_weight(const char *symbol) {
    for (int i = 0; i < PERIODIC_TABLE_LEN; i++) {
        if (strcmp(PERIODIC_TABLE[i].symbol, symbol) == 0)
            return PERIODIC_TABLE[i].weight;
    }
    return -1.0;
}

/* -------------------------------------------------------------------- */
/* Formula parsing                                                      */
/*                                                                      */
/* Case-insensitive with backtracking.  Element symbols in real         */
/* chemistry are an uppercase letter optionally followed by a lowercase  */
/* one (H, Cl, Na).  Users, however, often type sloppily ("c6h6",       */
/* "ccl4", "2nahco3").  To make that work we try to parse each position  */
/* as a 1- or 2-letter element and backtrack until the *whole* string    */
/* parses.                                                              */
/*                                                                      */
/* Preference order at each position:                                    */
/*   - If the input looks explicitly cased (uppercase then lowercase,    */
/*     e.g. "Co", "Cl") we prefer the 2-letter symbol, honouring intent. */
/*   - Otherwise (casual / lowercase) we prefer the 1-letter symbol and   */
/*     only fall back to the 2-letter one when needed.  This makes        */
/*     "nahco3" parse as Na-H-C-O3 rather than Na-H-Co3.                 */
/* -------------------------------------------------------------------- */
static int parse_number(const char *s, int *pos) {
    if (!isdigit((unsigned char)s[*pos])) return 1; /* default multiplier */
    int val = 0;
    while (isdigit((unsigned char)s[*pos])) {
        val = val * 10 + (s[*pos] - '0');
        (*pos)++;
    }
    return val;
}

/* Append text to the optional canonical-formula buffer (NULL = disabled). */
static void canon_append(char *canon, int *clen, const char *str) {
    if (!canon) return;
    while (*str && *clen < 126) canon[(*clen)++] = *str++;
    canon[*clen] = '\0';
}

/* Parse s[start..] up to end-of-string or a closing ')'.  On success
 * returns 1 with *weight set and *stop = index of the '\0' or ')'.
 * *far / far_sym track the deepest failure to build a helpful message.
 * If canon != NULL, the canonically-cased formula is written to it. */
static int parse_from(const char *s, int start, double *weight, int *stop,
                      int *far, char *far_sym, char *canon, int *clen) {
    if (s[start] == '\0' || s[start] == ')') {
        *weight = 0.0;
        *stop = start;
        return 1;
    }

    if (s[start] == '(') {
        double sub, rest;
        int inner_stop, rest_stop;
        int save = canon ? *clen : 0;
        canon_append(canon, clen, "(");
        if (!parse_from(s, start + 1, &sub, &inner_stop, far, far_sym, canon, clen)) {
            if (canon) { *clen = save; canon[*clen] = '\0'; }
            return 0;
        }
        if (s[inner_stop] != ')') {
            if (canon) { *clen = save; canon[*clen] = '\0'; }
            if (start > *far) { *far = start; far_sym[0] = '\0'; }
            return 0;
        }
        int p = inner_stop + 1;
        int mult = parse_number(s, &p);
        canon_append(canon, clen, ")");
        if (mult != 1) {
            char num[16]; snprintf(num, sizeof num, "%d", mult);
            canon_append(canon, clen, num);
        }
        if (!parse_from(s, p, &rest, &rest_stop, far, far_sym, canon, clen)) {
            if (canon) { *clen = save; canon[*clen] = '\0'; }
            return 0;
        }
        *weight = sub * mult + rest;
        *stop = rest_stop;
        return 1;
    }

    if (isalpha((unsigned char)s[start])) {
        char c0 = (char)toupper((unsigned char)s[start]);
        int next_lower = islower((unsigned char)s[start + 1]);
        int upper_first = isupper((unsigned char)s[start]);

        char one[3] = { c0, 0, 0 };
        char two[3] = { c0, next_lower ? (char)tolower((unsigned char)s[start + 1]) : 0, 0 };

        char cand[2][3];
        int  cand_len[2];
        int  ncand = 0;
        if (next_lower) {
            if (upper_first) {
                /* explicit "Co"/"Cl": prefer the 2-letter symbol */
                strcpy(cand[ncand], two); cand_len[ncand++] = 2;
                strcpy(cand[ncand], one); cand_len[ncand++] = 1;
            } else {
                /* casual / lowercase: prefer 1-letter, fall back to 2-letter */
                strcpy(cand[ncand], one); cand_len[ncand++] = 1;
                strcpy(cand[ncand], two); cand_len[ncand++] = 2;
            }
        } else {
            strcpy(cand[ncand], one); cand_len[ncand++] = 1;
        }

        for (int k = 0; k < ncand; k++) {
            double w = lookup_weight(cand[k]);
            if (w < 0) continue;
            int p = start + cand_len[k];
            int mult = parse_number(s, &p);
            double rest;
            int rest_stop;
            int save = canon ? *clen : 0;
            canon_append(canon, clen, cand[k]);
            if (mult != 1) {
                char num[16]; snprintf(num, sizeof num, "%d", mult);
                canon_append(canon, clen, num);
            }
            if (parse_from(s, p, &rest, &rest_stop, far, far_sym, canon, clen)) {
                *weight = w * mult + rest;
                *stop = rest_stop;
                return 1;
            }
            if (canon) { *clen = save; canon[*clen] = '\0'; }
        }

        /* Nothing worked here: remember it as an unknown-symbol guess. */
        if (start > *far) {
            *far = start;
            if (upper_first && next_lower) strcpy(far_sym, two);
            else strcpy(far_sym, one);
        }
        return 0;
    }

    if (start > *far) { *far = start; far_sym[0] = '\0'; }
    return 0;
}

/* Returns molecular weight, or -1.0 on error (bad_symbol filled if the
 * error was an unknown element).  If canonical != NULL it receives the
 * formula re-cased with canonical element symbols (e.g. "nahco3" ->
 * "NaHCO3"); the buffer must hold at least 128 bytes. */
static double molecular_weight(const char *formula, int *error,
                               char *bad_symbol, char *canonical) {
    *error = 0;
    if (bad_symbol) bad_symbol[0] = '\0';
    if (canonical) canonical[0] = '\0';

    if (formula[0] == '\0') { *error = 1; return -1.0; } /* empty formula */

    double w;
    int stop, far = -1, clen = 0;
    char far_sym[8] = {0};
    if (!parse_from(formula, 0, &w, &stop, &far, far_sym, canonical, &clen)) {
        *error = 1;
        if (bad_symbol) { strncpy(bad_symbol, far_sym, 7); bad_symbol[7] = '\0'; }
        return -1.0;
    }
    if (formula[stop] != '\0') { *error = 1; return -1.0; } /* stray ')' etc. */
    return w;
}

/* Split "2NaHCO3" / "1.5 NaOH" / "C6H6" into (equivalents, formula). */
static void split_equiv_formula(const char *in, double *equiv, char *formula, size_t fsz) {
    char *end;
    double e = strtod(in, &end);
    if (end != in && e > 0.0) {
        *equiv = e;
        while (*end == ' ') end++;
        strncpy(formula, end, fsz - 1);
    } else {
        *equiv = 1.0;
        strncpy(formula, in, fsz - 1);
    }
    formula[fsz - 1] = '\0';
}

/* -------------------------------------------------------------------- */
/* Colour                                                                */
/*                                                                      */
/* The `react` prompts run through three sections - reactants, then     */
/* products, then the name - and colour is what makes the switch        */
/* obvious as they scroll past.  One colour per section, switched off   */
/* unless stdout is a terminal, so redirected output stays plain.       */
/* NO_COLOR (any value) and TERM=dumb turn it off too.                  */
/* -------------------------------------------------------------------- */
static const char *col_react = "";   /* reactant prompts   */
static const char *col_prod  = "";   /* product prompts    */
static const char *col_id    = "";   /* reaction-name prompt */
static const char *col_dim   = "";   /* echoed results     */
static const char *col_off   = "";

static void init_colors(void) {
    const char *no_color = getenv("NO_COLOR");
    const char *term = getenv("TERM");
    if (no_color && no_color[0]) return;
    if (term && strcmp(term, "dumb") == 0) return;
    if (!isatty(STDOUT_FILENO)) return;

    col_react = "\033[1;36m";   /* cyan    */
    col_prod  = "\033[1;32m";   /* green   */
    col_id    = "\033[1;33m";   /* yellow  */
    col_dim   = "\033[2m";
    col_off   = "\033[0m";
}

/* -------------------------------------------------------------------- */
/* Units                                                                */
/* -------------------------------------------------------------------- */
/* Internally every mass is milligrams and every amount is millimoles:
 * that is the working scale of bench organic chemistry, so both the
 * stored records and the printed tables use it. */
static double unit_to_mg(double value, const char *unit) {
    if (strcmp(unit, "mg") == 0) return value;
    if (strcmp(unit, "g") == 0) return value * 1e3;
    if (strcmp(unit, "kg") == 0) return value * 1e6;
    if (strcmp(unit, "ug") == 0) return value * 1e-3;
    return -1.0; /* unsupported */
}

/* -------------------------------------------------------------------- */
/* Data model                                                           */
/* -------------------------------------------------------------------- */
#define MAX_REACTANTS 10
#define MAX_PRODUCTS  10
#define MAX_REACTIONS 1000
#define MAX_ID        32   /* ids may be user-chosen names, not just hex */

typedef struct {
    char label;
    char formula[64];
    double equiv;      /* stoichiometric coefficient (equivalents) */
    double mass_mg;
    char unit[8];      /* unit the mass was typed in (kept for the record) */
    double mw;
    double mmol;
} Reactant;

typedef struct {
    char formula[64];
    double equiv;      /* stoichiometric coefficient (equivalents) */
    double mw;
    double theoretical_mmol;
    double theoretical_yield_mg;
    int has_actual;
    double actual_yield_mg;
    double percent_yield;
} Product;

typedef struct {
    char id[MAX_ID + 1];
    char created[32];
    Reactant reactants[MAX_REACTANTS];
    int nreactants;
    Product products[MAX_PRODUCTS];
    int nproducts;
    char limiting_label;
} Reaction;

/* -------------------------------------------------------------------- */
/* Storage: ~/.local/share/chemcalc/reactions.dat                       */
/*                                                                      */
/* First line is a version marker; then one line per reaction, with     */
/* tab-delimited fields:                                                */
/*   id  created  limiting  reactants  products                         */
/* The id is a random hex string or a user-chosen name; valid_id keeps   */
/* it clear of the separators used here.                                 */
/* reactants: ';'-joined "label,formula,equiv,mass_mg,unit,mw,mmol"     */
/* products:  ';'-joined                                                */
/*   "formula,equiv,mw,theo_mmol,theo_mg,has_actual,actual_mg,percent"  */
/*                                                                      */
/* v1 files carry no marker and hold grams/moles instead of mg/mmol;    */
/* they are converted while loading and rewritten as v2 on the next     */
/* save.  The reactant 'unit' column only records what the user typed;  */
/* the mass itself is always milligrams (in v1 it was in that unit).    */
/* -------------------------------------------------------------------- */
#define STORE_VERSION_LINE "#chemcalc v2"
static void store_path(char *out, size_t outsz) {
    const char *home = getenv("HOME");
    if (!home) home = ".";
    snprintf(out, outsz, "%s/.local/share/chemcalc", home);
}

static void ensure_store_dir(void) {
    char dir[512], tmp[512];
    store_path(dir, sizeof(dir));
    snprintf(tmp, sizeof(tmp), "%s", dir);
    for (char *p = tmp + 1; *p; p++) {
        if (*p == '/') {
            *p = '\0';
            mkdir(tmp, 0755);
            *p = '/';
        }
    }
    mkdir(tmp, 0755);
}

static void store_file(char *out, size_t outsz) {
    const char *home = getenv("HOME");
    if (!home) home = ".";
    snprintf(out, outsz, "%s/.local/share/chemcalc/reactions.dat", home);
}

static int load_reactions(Reaction *reactions, int max) {
    char path[512];
    store_file(path, sizeof(path));
    FILE *f = fopen(path, "r");
    if (!f) return 0;

    char line[8192];
    int n = 0;

    /* A missing marker means a v1 (grams/moles) file: rewind and convert. */
    int version = 1;
    if (fgets(line, sizeof(line), f)) {
        line[strcspn(line, "\n")] = '\0';
        if (strcmp(line, STORE_VERSION_LINE) == 0) version = 2;
    }
    if (version == 1) rewind(f);

    while (fgets(line, sizeof(line), f) && n < max) {
        line[strcspn(line, "\n")] = '\0';
        if (line[0] == '\0') continue;

        Reaction *r = &reactions[n];
        memset(r, 0, sizeof(Reaction));

        char *saveptr;
        char *tok;
        int field = 0;
        char reactants_field[4096] = "";
        char products_field[4096] = "";

        tok = strtok_r(line, "\t", &saveptr);
        while (tok) {
            switch (field) {
                case 0: strncpy(r->id, tok, sizeof(r->id) - 1); break;
                case 1: strncpy(r->created, tok, sizeof(r->created) - 1); break;
                case 2: r->limiting_label = tok[0]; break;
                case 3: strncpy(reactants_field, tok, sizeof(reactants_field) - 1); break;
                case 4: strncpy(products_field, tok, sizeof(products_field) - 1); break;
            }
            field++;
            tok = strtok_r(NULL, "\t", &saveptr);
        }

        /* reactants */
        r->nreactants = 0;
        char *rsave;
        char *rec = strtok_r(reactants_field, ";", &rsave);
        while (rec && r->nreactants < MAX_REACTANTS) {
            Reactant *rt = &r->reactants[r->nreactants];
            char buf[512];
            strncpy(buf, rec, sizeof(buf) - 1);
            buf[sizeof(buf) - 1] = '\0';

            char *csave;
            char *c = strtok_r(buf, ",", &csave);
            int cf = 0;
            while (c) {
                switch (cf) {
                    case 0: rt->label = c[0]; break;
                    case 1: strncpy(rt->formula, c, sizeof(rt->formula) - 1); break;
                    case 2: rt->equiv = atof(c); break;
                    case 3: rt->mass_mg = atof(c); break;
                    case 4: strncpy(rt->unit, c, sizeof(rt->unit) - 1); break;
                    case 5: rt->mw = atof(c); break;
                    case 6: rt->mmol = atof(c); break;
                }
                cf++;
                c = strtok_r(NULL, ",", &csave);
            }
            if (rt->equiv <= 0.0) rt->equiv = 1.0;
            if (version == 1) {
                /* v1 held the mass in rt->unit and the amount in moles. */
                double mg = unit_to_mg(rt->mass_mg, rt->unit);
                if (mg >= 0.0) rt->mass_mg = mg;
                rt->mmol *= 1e3;
            }
            r->nreactants++;
            rec = strtok_r(NULL, ";", &rsave);
        }

        /* products */
        r->nproducts = 0;
        char *psave;
        char *prec = strtok_r(products_field, ";", &psave);
        while (prec && r->nproducts < MAX_PRODUCTS) {
            Product *p = &r->products[r->nproducts];
            char buf[512];
            strncpy(buf, prec, sizeof(buf) - 1);
            buf[sizeof(buf) - 1] = '\0';

            char *csave;
            char *c = strtok_r(buf, ",", &csave);
            int cf = 0;
            while (c) {
                switch (cf) {
                    case 0: strncpy(p->formula, c, sizeof(p->formula) - 1); break;
                    case 1: p->equiv = atof(c); break;
                    case 2: p->mw = atof(c); break;
                    case 3: p->theoretical_mmol = atof(c); break;
                    case 4: p->theoretical_yield_mg = atof(c); break;
                    case 5: p->has_actual = atoi(c); break;
                    case 6: p->actual_yield_mg = atof(c); break;
                    case 7: p->percent_yield = atof(c); break;
                }
                cf++;
                c = strtok_r(NULL, ",", &csave);
            }
            if (p->equiv <= 0.0) p->equiv = 1.0;
            if (version == 1) {
                /* v1 held moles and grams. */
                p->theoretical_mmol *= 1e3;
                p->theoretical_yield_mg *= 1e3;
                p->actual_yield_mg *= 1e3;
            }
            r->nproducts++;
            prec = strtok_r(NULL, ";", &psave);
        }

        n++;
    }
    fclose(f);
    return n;
}

static void save_reactions(Reaction *reactions, int n) {
    ensure_store_dir();
    char path[512];
    store_file(path, sizeof(path));
    FILE *f = fopen(path, "w");
    if (!f) {
        fprintf(stderr, "Error: could not write to %s\n", path);
        exit(1);
    }
    fprintf(f, "%s\n", STORE_VERSION_LINE);
    for (int i = 0; i < n; i++) {
        Reaction *r = &reactions[i];
        fprintf(f, "%s\t%s\t%c\t",
                r->id, r->created, r->limiting_label ? r->limiting_label : '-');
        for (int j = 0; j < r->nreactants; j++) {
            Reactant *rt = &r->reactants[j];
            fprintf(f, "%c,%s,%.6f,%.6f,%s,%.6f,%.6f",
                    rt->label, rt->formula, rt->equiv, rt->mass_mg, rt->unit,
                    rt->mw, rt->mmol);
            if (j != r->nreactants - 1) fprintf(f, ";");
        }
        fprintf(f, "\t");
        for (int j = 0; j < r->nproducts; j++) {
            Product *p = &r->products[j];
            fprintf(f, "%s,%.6f,%.6f,%.6f,%.6f,%d,%.6f,%.6f",
                    p->formula, p->equiv, p->mw, p->theoretical_mmol,
                    p->theoretical_yield_mg, p->has_actual,
                    p->actual_yield_mg, p->percent_yield);
            if (j != r->nproducts - 1) fprintf(f, ";");
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

/* Look up a reaction by ID.  An exact ID always wins; otherwise a unique
 * prefix match is accepted.  When a prefix matches more than one reaction
 * *ambiguous is set (if provided) and NULL is returned. */
static Reaction *find_reaction(Reaction *reactions, int n, const char *id,
                               int *ambiguous) {
    if (ambiguous) *ambiguous = 0;

    for (int i = 0; i < n; i++)
        if (strcmp(reactions[i].id, id) == 0) return &reactions[i];

    Reaction *match = NULL;
    int count = 0;
    size_t len = strlen(id);
    for (int i = 0; i < n; i++) {
        if (strncmp(reactions[i].id, id, len) == 0) { match = &reactions[i]; count++; }
    }
    if (count == 1) return match;
    if (count > 1 && ambiguous) *ambiguous = 1;
    return NULL;
}

static int id_taken(Reaction *reactions, int n, const char *id) {
    for (int i = 0; i < n; i++)
        if (strcmp(reactions[i].id, id) == 0) return 1;
    return 0;
}

/* Ids double as names, so they are restricted to characters that stay easy
 * to type on the command line and can never collide with the store's field
 * separators (tab, ';' and ','). */
static int valid_id(const char *id, char *err, size_t errsz) {
    size_t len = strlen(id);
    if (len == 0) {
        snprintf(err, errsz, "name must not be empty");
        return 0;
    }
    if (len > MAX_ID) {
        snprintf(err, errsz, "name must be at most %d characters", MAX_ID);
        return 0;
    }
    for (const char *p = id; *p; p++)
        if (!isalnum((unsigned char)*p) && *p != '-' && *p != '_' && *p != '.') {
            snprintf(err, errsz,
                     "name may only contain letters, digits, '-', '_' and '.'");
            return 0;
        }
    return 1;
}

static void gen_id(Reaction *reactions, int n, char *out) {
    static int seeded = 0;
    if (!seeded) { srand((unsigned)(time(NULL) ^ getpid())); seeded = 1; }
    const char *hex = "0123456789abcdef";
    do {
        for (int i = 0; i < 8; i++) out[i] = hex[rand() % 16];
        out[8] = '\0';
    } while (id_taken(reactions, n, out));
}

static void now_iso(char *out, size_t outsz) {
    time_t t = time(NULL);
    struct tm tmv;
    localtime_r(&t, &tmv);
    strftime(out, outsz, "%Y-%m-%dT%H:%M:%S", &tmv);
}

/* -------------------------------------------------------------------- */
/* Simple table printer                                                 */
/* -------------------------------------------------------------------- */
#define MAX_COLS 24
#define MAX_ROWS 256

static void print_table(const char *headers[], int ncols,
                         char rows[][MAX_COLS][128], int nrows) {
    int widths[MAX_COLS];
    for (int c = 0; c < ncols; c++) widths[c] = (int)strlen(headers[c]);
    for (int r = 0; r < nrows; r++)
        for (int c = 0; c < ncols; c++) {
            int len = (int)strlen(rows[r][c]);
            if (len > widths[c]) widths[c] = len;
        }

    void print_sep(void) {
        for (int c = 0; c < ncols; c++) {
            putchar('+');
            for (int i = 0; i < widths[c] + 2; i++) putchar('-');
        }
        printf("+\n");
    }

    print_sep();
    printf("|");
    for (int c = 0; c < ncols; c++) printf(" %-*s |", widths[c], headers[c]);
    printf("\n");
    print_sep();
    for (int r = 0; r < nrows; r++) {
        printf("|");
        for (int c = 0; c < ncols; c++) printf(" %-*s |", widths[c], rows[r][c]);
        printf("\n");
    }
    print_sep();
}

/* -------------------------------------------------------------------- */
/* Helpers for displaying a reaction                                     */
/* -------------------------------------------------------------------- */
static void species_name(char *out, size_t sz, double equiv, const char *formula) {
    if (equiv == 1.0) snprintf(out, sz, "%s", formula);
    else              snprintf(out, sz, "%g %s", equiv, formula);
}

static void build_equation(Reaction *r, char *out, size_t sz) {
    out[0] = '\0';
    for (int i = 0; i < r->nreactants; i++) {
        char nm[96];
        species_name(nm, sizeof nm, r->reactants[i].equiv, r->reactants[i].formula);
        strncat(out, nm, sz - strlen(out) - 1);
        if (i != r->nreactants - 1) strncat(out, " + ", sz - strlen(out) - 1);
    }
    strncat(out, " -> ", sz - strlen(out) - 1);
    for (int i = 0; i < r->nproducts; i++) {
        char nm[96];
        species_name(nm, sizeof nm, r->products[i].equiv, r->products[i].formula);
        strncat(out, nm, sz - strlen(out) - 1);
        if (i != r->nproducts - 1) strncat(out, " + ", sz - strlen(out) - 1);
    }
}

/* -------------------------------------------------------------------- */
/* Printing a single reaction (used by react/show/yield)                */
/*                                                                      */
/* The table is transposed: the header row holds the species names,     */
/* then one row each for molecular weight, mass, amount, the actual     */
/* amount obtained (mg and mmol) and % yield.  Masses are milligrams    */
/* and amounts millimoles throughout.                                   */
/* -------------------------------------------------------------------- */
static void print_reaction(Reaction *r) {
    int nspec = r->nreactants + r->nproducts;
    int ncols = 1 + nspec;
    if (ncols > MAX_COLS) ncols = MAX_COLS;

    char header_cells[MAX_COLS][128];
    const char *headers[MAX_COLS];
    snprintf(header_cells[0], 128, "%s", "");

    int c = 1;
    for (int i = 0; i < r->nreactants && c < MAX_COLS; i++, c++) {
        char nm[96];
        species_name(nm, sizeof nm, r->reactants[i].equiv, r->reactants[i].formula);
        snprintf(header_cells[c], 128, "%c:%s%s", r->reactants[i].label, nm,
                 r->reactants[i].label == r->limiting_label ? " *" : "");
    }
    for (int i = 0; i < r->nproducts && c < MAX_COLS; i++, c++) {
        char nm[96];
        species_name(nm, sizeof nm, r->products[i].equiv, r->products[i].formula);
        snprintf(header_cells[c], 128, "P%d:%s", i + 1, nm);
    }
    for (int k = 0; k < ncols; k++) headers[k] = header_cells[k];

    char rows[8][MAX_COLS][128];
    int nrows = 0;

    /* Row: molecular weight */
    snprintf(rows[nrows][0], 128, "MW (g/mol)");
    c = 1;
    for (int i = 0; i < r->nreactants && c < MAX_COLS; i++, c++)
        snprintf(rows[nrows][c], 128, "%.3f", r->reactants[i].mw);
    for (int i = 0; i < r->nproducts && c < MAX_COLS; i++, c++)
        snprintf(rows[nrows][c], 128, "%.3f", r->products[i].mw);
    nrows++;

    /* Row: mass (mg) */
    snprintf(rows[nrows][0], 128, "Mass (mg)");
    c = 1;
    for (int i = 0; i < r->nreactants && c < MAX_COLS; i++, c++)
        snprintf(rows[nrows][c], 128, "%.3f", r->reactants[i].mass_mg);
    for (int i = 0; i < r->nproducts && c < MAX_COLS; i++, c++)
        snprintf(rows[nrows][c], 128, "%.3f", r->products[i].theoretical_yield_mg);
    nrows++;

    /* Row: amount (mmol) */
    snprintf(rows[nrows][0], 128, "Moles (mmol)");
    c = 1;
    for (int i = 0; i < r->nreactants && c < MAX_COLS; i++, c++)
        snprintf(rows[nrows][c], 128, "%.4f", r->reactants[i].mmol);
    for (int i = 0; i < r->nproducts && c < MAX_COLS; i++, c++)
        snprintf(rows[nrows][c], 128, "%.4f", r->products[i].theoretical_mmol);
    nrows++;

    /* Row: mmol relative to reactant A, the form the amounts are entered in */
    snprintf(rows[nrows][0], 128, "Equiv (vs %c)",
             r->nreactants ? r->reactants[0].label : 'A');
    c = 1;
    double mmol_a = r->nreactants ? r->reactants[0].mmol : 0.0;
    for (int i = 0; i < r->nreactants && c < MAX_COLS; i++, c++) {
        if (mmol_a > 0.0) snprintf(rows[nrows][c], 128, "%.4g",
                                   r->reactants[i].mmol / mmol_a);
        else              snprintf(rows[nrows][c], 128, "-");
    }
    for (int i = 0; i < r->nproducts && c < MAX_COLS; i++, c++)
        snprintf(rows[nrows][c], 128, "-");
    nrows++;

    /* Rows: the amount actually obtained, as recorded by the yield command */
    snprintf(rows[nrows][0], 128, "Actual (mg)");
    c = 1;
    for (int i = 0; i < r->nreactants && c < MAX_COLS; i++, c++)
        snprintf(rows[nrows][c], 128, "-");
    for (int i = 0; i < r->nproducts && c < MAX_COLS; i++, c++) {
        if (r->products[i].has_actual)
            snprintf(rows[nrows][c], 128, "%.3f", r->products[i].actual_yield_mg);
        else
            snprintf(rows[nrows][c], 128, "-");
    }
    nrows++;

    snprintf(rows[nrows][0], 128, "Actual (mmol)");
    c = 1;
    for (int i = 0; i < r->nreactants && c < MAX_COLS; i++, c++)
        snprintf(rows[nrows][c], 128, "-");
    for (int i = 0; i < r->nproducts && c < MAX_COLS; i++, c++) {
        Product *p = &r->products[i];
        if (p->has_actual && p->mw > 0.0)
            snprintf(rows[nrows][c], 128, "%.4f", p->actual_yield_mg / p->mw);
        else
            snprintf(rows[nrows][c], 128, "-");
    }
    nrows++;

    /* Row: yield (%) */
    snprintf(rows[nrows][0], 128, "Yield (%%)");
    c = 1;
    for (int i = 0; i < r->nreactants && c < MAX_COLS; i++, c++)
        snprintf(rows[nrows][c], 128, "-");
    for (int i = 0; i < r->nproducts && c < MAX_COLS; i++, c++) {
        if (r->products[i].has_actual)
            snprintf(rows[nrows][c], 128, "%.2f", r->products[i].percent_yield);
        else
            snprintf(rows[nrows][c], 128, "-");
    }
    nrows++;

    char eq[512];
    build_equation(r, eq, sizeof eq);
    printf("Equation: %s\n", eq);
    if (r->limiting_label) {
        for (int i = 0; i < r->nreactants; i++)
            if (r->reactants[i].label == r->limiting_label)
                printf("Limiting reagent: %c (%s)\n", r->limiting_label,
                       r->reactants[i].formula);
    }
    printf("\n");

    print_table(headers, ncols, rows, nrows);
    printf("\n(Product Mass/Moles are theoretical; Yield %% needs a recorded actual.\n");
    printf(" '*' marks the limiting reagent.)\n");
}

/* -------------------------------------------------------------------- */
/* Commands                                                              */
/* -------------------------------------------------------------------- */
static void cmd_mw(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: chemcalc mw FORMULA [--mass VALUE] [--unit g|mg|kg|ug]\n");
        exit(1);
    }
    /* A leading coefficient (e.g. "2H2O") is accepted but ignored here:
       molecular weight is a per-species quantity. */
    double equiv;
    char formula[128];
    split_equiv_formula(argv[2], &equiv, formula, sizeof(formula));

    double mass = -1.0;
    const char *unit = "mg";

    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--mass") == 0 && i + 1 < argc) mass = atof(argv[++i]);
        else if (strcmp(argv[i], "--unit") == 0 && i + 1 < argc) unit = argv[++i];
    }

    int error = 0;
    char bad_symbol[8], canonical[128];
    double mw = molecular_weight(formula, &error, bad_symbol, canonical);
    if (error) {
        if (formula[0] == '\0') fprintf(stderr, "Error: no formula given\n");
        else if (bad_symbol[0]) fprintf(stderr, "Error: unknown element '%s'\n", bad_symbol);
        else fprintf(stderr, "Error: could not parse formula '%s'\n", formula);
        exit(1);
    }

    printf("Formula: %s\n", canonical);
    printf("Molecular weight: %.4f g/mol\n", mw);

    if (mass >= 0.0) {
        double mass_mg = unit_to_mg(mass, unit);
        if (mass_mg < 0.0) {
            fprintf(stderr, "Error: unsupported unit '%s' (use g/mg/kg/ug)\n", unit);
            exit(1);
        }
        double mmol = mass_mg / mw;
        printf("Mass: %g %s -> %.4f mmol\n", mass, unit, mmol);
    }
}

static void trim_newline(char *s) { s[strcspn(s, "\n")] = '\0'; }

static void cmd_react(void) {
    Reaction r;
    memset(&r, 0, sizeof(r));
    char line[256];

    printf("%s== Reactants ==%s Blank formula to finish (need at least 1).\n",
           col_react, col_off);
    printf("Prefix a coefficient for equivalents, e.g. 2NaHCO3. Case is ignored.\n");
    printf("Reactant A is weighed out; the rest are given as a multiple of A's mmol.\n");

    char label = 'A';
    while (r.nreactants < MAX_REACTANTS) {
        printf("%sReactant %c formula%s (e.g. 2NaHCO3), or blank to stop: ",
               col_react, label, col_off);
        fflush(stdout);
        if (!fgets(line, sizeof(line), stdin)) break;
        trim_newline(line);
        if (line[0] == '\0') break;

        Reactant *rt = &r.reactants[r.nreactants];
        rt->label = label;
        split_equiv_formula(line, &rt->equiv, rt->formula, sizeof(rt->formula));

        int error = 0;
        char bad_symbol[8], canonical[128];
        rt->mw = molecular_weight(rt->formula, &error, bad_symbol, canonical);
        if (error) {
            if (bad_symbol[0]) fprintf(stderr, "Error: unknown element '%s'\n", bad_symbol);
            else fprintf(stderr, "Error: could not parse formula '%s'\n", rt->formula);
            exit(1);
        }
        snprintf(rt->formula, sizeof(rt->formula), "%s", canonical);

        if (r.nreactants == 0) {
            /* Reactant A sets the scale of the run: it is weighed out. */
            printf("%sReactant %c mass value%s: ", col_react, label, col_off);
            fflush(stdout);
            if (!fgets(line, sizeof(line), stdin)) exit(1);
            trim_newline(line);
            double entered = atof(line);

            printf("%sReactant %c unit%s [g/mg] (default mg): ",
                   col_react, label, col_off);
            fflush(stdout);
            if (!fgets(line, sizeof(line), stdin)) exit(1);
            trim_newline(line);
            snprintf(rt->unit, sizeof(rt->unit), "%.7s", line[0] ? line : "mg");

            rt->mass_mg = unit_to_mg(entered, rt->unit);
            if (rt->mass_mg < 0.0) {
                fprintf(stderr, "Error: unsupported unit '%s'\n", rt->unit);
                exit(1);
            }
            rt->mmol = rt->mass_mg / rt->mw;
        } else {
            /* Everything after A is entered as "how many times A's mmol",
             * and the mass needed to hit that amount is worked out here.
             * The default is the ratio the coefficients already imply, so
             * C16H10 + 4Br2 gives B 4x the mmol of A. A bare number is a
             * ratio; a number with a unit ("500 mg") is a mass instead. */
            double mmol_a = r.reactants[0].mmol;
            double def_ratio = r.reactants[0].equiv > 0.0
                             ? rt->equiv / r.reactants[0].equiv : rt->equiv;

            printf("%sReactant %c mmol relative to %c%s [default %g x], "
                   "or a mass like \"500 mg\": ",
                   col_react, label, r.reactants[0].label, col_off, def_ratio);
            fflush(stdout);
            if (!fgets(line, sizeof(line), stdin)) exit(1);
            trim_newline(line);

            char *end;
            double val = strtod(line, &end);
            int have_number = (end != line);
            while (*end == ' ' || *end == '\t') end++;

            if (line[0] == '\0' || (have_number && *end == '\0')) {
                /* Ratio against A (blank keeps the stoichiometric one). */
                double ratio = line[0] == '\0' ? def_ratio : val;
                if (ratio <= 0.0) {
                    fprintf(stderr, "Error: ratio must be greater than 0\n");
                    exit(1);
                }
                rt->mmol = mmol_a * ratio;
                rt->mass_mg = rt->mmol * rt->mw;
                snprintf(rt->unit, sizeof(rt->unit), "mg");
                printf("%s  -> %.4f mmol, %.3f mg%s\n",
                       col_dim, rt->mmol, rt->mass_mg, col_off);
            } else if (have_number) {
                /* An explicit mass, for a reagent that was weighed out. */
                snprintf(rt->unit, sizeof(rt->unit), "%.7s", end);
                rt->mass_mg = unit_to_mg(val, rt->unit);
                if (rt->mass_mg < 0.0) {
                    fprintf(stderr, "Error: unsupported unit '%s'\n", rt->unit);
                    exit(1);
                }
                rt->mmol = rt->mass_mg / rt->mw;
                printf("%s  -> %.4f mmol (%.4g x %c)%s\n", col_dim, rt->mmol,
                       mmol_a > 0.0 ? rt->mmol / mmol_a : 0.0,
                       r.reactants[0].label, col_off);
            } else {
                fprintf(stderr, "Error: could not read '%s' as a ratio or a mass\n",
                        line);
                exit(1);
            }
        }

        r.nreactants++;
        label++;
    }

    if (r.nreactants == 0) {
        fprintf(stderr, "No reactants entered, aborting.\n");
        exit(1);
    }

    printf("\n%s== Products ==%s Blank formula to finish (need at least 1).\n",
           col_prod, col_off);
    while (r.nproducts < MAX_PRODUCTS) {
        printf("%sProduct %d formula%s (e.g. C6H6 or 2H2O), or blank to stop: ",
               col_prod, r.nproducts + 1, col_off);
        fflush(stdout);
        if (!fgets(line, sizeof(line), stdin)) break;
        trim_newline(line);
        if (line[0] == '\0') break;

        Product *p = &r.products[r.nproducts];
        split_equiv_formula(line, &p->equiv, p->formula, sizeof(p->formula));

        int error = 0;
        char bad_symbol[8], canonical[128];
        p->mw = molecular_weight(p->formula, &error, bad_symbol, canonical);
        if (error) {
            if (bad_symbol[0]) fprintf(stderr, "Error: unknown element '%s'\n", bad_symbol);
            else fprintf(stderr, "Error: could not parse formula '%s'\n", p->formula);
            exit(1);
        }
        snprintf(p->formula, sizeof(p->formula), "%s", canonical);
        r.nproducts++;
    }

    if (r.nproducts == 0) {
        fprintf(stderr, "No products entered, aborting.\n");
        exit(1);
    }

    /* Limiting reagent = smallest mmol/equiv; the reaction extent is that
     * ratio, and each product's theoretical amount scales by its equiv. */
    int lim = 0;
    for (int i = 1; i < r.nreactants; i++)
        if (r.reactants[i].mmol / r.reactants[i].equiv <
            r.reactants[lim].mmol / r.reactants[lim].equiv) lim = i;
    r.limiting_label = r.reactants[lim].label;
    double extent = r.reactants[lim].mmol / r.reactants[lim].equiv;

    for (int i = 0; i < r.nproducts; i++) {
        r.products[i].theoretical_mmol = extent * r.products[i].equiv;
        r.products[i].theoretical_yield_mg =
            r.products[i].theoretical_mmol * r.products[i].mw;
        r.products[i].has_actual = 0;
    }

    Reaction all[MAX_REACTIONS];
    int n = load_reactions(all, MAX_REACTIONS);
    if (n >= MAX_REACTIONS) {
        fprintf(stderr, "Error: reaction store full\n");
        exit(1);
    }

    /* The id may be a name of the user's choosing; blank picks a random one. */
    printf("\n%s== Name ==%s The id you will look this reaction up by.\n",
           col_id, col_off);
    for (;;) {
        printf("%sReaction name%s (blank for a random id): ", col_id, col_off);
        fflush(stdout);
        if (!fgets(line, sizeof(line), stdin)) { gen_id(all, n, r.id); break; }
        trim_newline(line);
        if (line[0] == '\0') { gen_id(all, n, r.id); break; }

        char err[128];
        if (!valid_id(line, err, sizeof err)) {
            fprintf(stderr, "Error: %s\n", err);
            continue;
        }
        if (id_taken(all, n, line)) {
            fprintf(stderr, "Error: id '%s' is already in use\n", line);
            continue;
        }
        snprintf(r.id, sizeof(r.id), "%.*s", MAX_ID, line);
        break;
    }
    now_iso(r.created, sizeof(r.created));

    all[n++] = r;
    save_reactions(all, n);

    print_reaction(&r);
    printf("\nSaved as reaction id: %s%s%s\n", col_id, r.id, col_off);
}

static void cmd_list(void) {
    Reaction all[MAX_REACTIONS];
    int n = load_reactions(all, MAX_REACTIONS);
    if (n == 0) {
        printf("No saved reactions yet.\n");
        return;
    }

    const char *headers[] = {"ID", "Created", "Reactants", "Products",
                             "Theo. yield (mg)", "% yield"};
    char rows[MAX_ROWS][MAX_COLS][128];
    int nrows = n < MAX_ROWS ? n : MAX_ROWS;

    for (int i = 0; i < nrows; i++) {
        Reaction *r = &all[i];
        char reactant_str[256] = "";
        for (int j = 0; j < r->nreactants; j++) {
            char nm[96];
            species_name(nm, sizeof nm, r->reactants[j].equiv, r->reactants[j].formula);
            strncat(reactant_str, nm, sizeof(reactant_str) - strlen(reactant_str) - 1);
            if (j != r->nreactants - 1)
                strncat(reactant_str, " + ", sizeof(reactant_str) - strlen(reactant_str) - 1);
        }
        char product_str[256] = "";
        double theo_total = 0.0;
        for (int j = 0; j < r->nproducts; j++) {
            char nm[96];
            species_name(nm, sizeof nm, r->products[j].equiv, r->products[j].formula);
            strncat(product_str, nm, sizeof(product_str) - strlen(product_str) - 1);
            if (j != r->nproducts - 1)
                strncat(product_str, " + ", sizeof(product_str) - strlen(product_str) - 1);
            theo_total += r->products[j].theoretical_yield_mg;
        }
        snprintf(rows[i][0], 128, "%.*s", MAX_ID, r->id);
        snprintf(rows[i][1], 128, "%s", r->created);
        snprintf(rows[i][2], 128, "%s", reactant_str);
        snprintf(rows[i][3], 128, "%s", product_str);
        snprintf(rows[i][4], 128, "%.3f", theo_total);
        if (r->nproducts > 0 && r->products[0].has_actual)
            snprintf(rows[i][5], 128, "%.2f", r->products[0].percent_yield);
        else
            snprintf(rows[i][5], 128, "-");
    }
    print_table(headers, 6, rows, nrows);
}

static void cmd_delete(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: chemcalc delete ID [-f|--force]\n");
        exit(1);
    }
    int force = 0;
    for (int i = 3; i < argc; i++)
        if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--force") == 0) force = 1;

    Reaction all[MAX_REACTIONS];
    int n = load_reactions(all, MAX_REACTIONS);
    int ambiguous = 0;
    Reaction *r = find_reaction(all, n, argv[2], &ambiguous);
    if (!r) {
        if (ambiguous)
            fprintf(stderr, "Ambiguous id '%s' matches more than one reaction\n", argv[2]);
        else
            fprintf(stderr, "No reaction found with id '%s'\n", argv[2]);
        exit(1);
    }

    char id[MAX_ID + 1], eq[512];
    snprintf(id, sizeof(id), "%s", r->id);
    build_equation(r, eq, sizeof eq);

    /* Deleting is not undoable, so confirm unless told not to.  Without a
     * terminal to ask on, --force is required. */
    if (!force) {
        if (!isatty(STDIN_FILENO)) {
            fprintf(stderr, "Refusing to delete '%s' without a terminal; "
                            "pass --force.\n", id);
            exit(1);
        }
        printf("Delete reaction %s (%s)? [y/N]: ", id, eq);
        fflush(stdout);
        char line[64];
        if (!fgets(line, sizeof(line), stdin)) { printf("\nCancelled.\n"); return; }
        trim_newline(line);
        if (line[0] != 'y' && line[0] != 'Y') {
            printf("Cancelled.\n");
            return;
        }
    }

    int idx = (int)(r - all);
    for (int i = idx; i < n - 1; i++) all[i] = all[i + 1];
    n--;
    save_reactions(all, n);
    printf("Deleted reaction %s (%s)\n", id, eq);
}

static void cmd_rename(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage: chemcalc rename ID NEW_NAME\n");
        exit(1);
    }
    const char *new_id = argv[3];

    char err[128];
    if (!valid_id(new_id, err, sizeof err)) {
        fprintf(stderr, "Error: %s\n", err);
        exit(1);
    }

    Reaction all[MAX_REACTIONS];
    int n = load_reactions(all, MAX_REACTIONS);
    int ambiguous = 0;
    Reaction *r = find_reaction(all, n, argv[2], &ambiguous);
    if (!r) {
        if (ambiguous)
            fprintf(stderr, "Ambiguous id '%s' matches more than one reaction\n", argv[2]);
        else
            fprintf(stderr, "No reaction found with id '%s'\n", argv[2]);
        exit(1);
    }
    if (strcmp(r->id, new_id) == 0) {
        printf("Reaction is already named %s\n", new_id);
        return;
    }
    if (id_taken(all, n, new_id)) {
        fprintf(stderr, "Error: id '%s' is already in use\n", new_id);
        exit(1);
    }

    char old_id[MAX_ID + 1];
    snprintf(old_id, sizeof(old_id), "%s", r->id);
    snprintf(r->id, sizeof(r->id), "%s", new_id);
    save_reactions(all, n);
    printf("Renamed %s -> %s\n", old_id, r->id);
}

static void cmd_show(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: chemcalc show ID\n");
        exit(1);
    }
    Reaction all[MAX_REACTIONS];
    int n = load_reactions(all, MAX_REACTIONS);
    int ambiguous = 0;
    Reaction *r = find_reaction(all, n, argv[2], &ambiguous);
    if (!r) {
        if (ambiguous)
            fprintf(stderr, "Ambiguous id '%s' matches more than one reaction\n", argv[2]);
        else
            fprintf(stderr, "No reaction found with id '%s'\n", argv[2]);
        exit(1);
    }
    print_reaction(r);
}

static void cmd_yield(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage: chemcalc yield ID ACTUAL_VALUE "
                        "[--unit g|mg|kg|ug] [--product N]\n");
        exit(1);
    }
    const char *id = argv[2];
    double actual = atof(argv[3]);
    const char *unit = "mg";
    int pidx = 1;
    for (int i = 4; i < argc; i++) {
        if (strcmp(argv[i], "--unit") == 0 && i + 1 < argc) unit = argv[++i];
        else if (strcmp(argv[i], "--product") == 0 && i + 1 < argc) pidx = atoi(argv[++i]);
    }

    Reaction all[MAX_REACTIONS];
    int n = load_reactions(all, MAX_REACTIONS);
    int ambiguous = 0;
    Reaction *r = find_reaction(all, n, id, &ambiguous);
    if (!r) {
        if (ambiguous)
            fprintf(stderr, "Ambiguous id '%s' matches more than one reaction\n", id);
        else
            fprintf(stderr, "No reaction found with id '%s'\n", id);
        exit(1);
    }
    if (pidx < 1 || pidx > r->nproducts) {
        fprintf(stderr, "Error: product %d out of range (1..%d)\n", pidx, r->nproducts);
        exit(1);
    }

    Product *p = &r->products[pidx - 1];
    double actual_mg = unit_to_mg(actual, unit);
    if (actual_mg < 0.0) {
        fprintf(stderr, "Error: unsupported unit '%s'\n", unit);
        exit(1);
    }
    p->has_actual = 1;
    p->actual_yield_mg = actual_mg;
    p->percent_yield = p->theoretical_yield_mg > 0
        ? (actual_mg / p->theoretical_yield_mg) * 100.0 : 0.0;

    save_reactions(all, n);
    print_reaction(r);
}

static void print_usage(void) {
    printf("chemcalc - simple chemistry calculator CLI\n\n");
    printf("Formulas are case-insensitive (C6H6, c6h6, ccl4).\n");
    printf("Reactants/products may carry a coefficient (equivalents), e.g. 2NaHCO3.\n");
    printf("Masses default to mg and amounts are reported in mmol.\n\n");
    printf("Usage:\n");
    printf("  chemcalc mw FORMULA [--mass VALUE] [--unit g|mg|kg|ug]\n");
    printf("  chemcalc react\n");
    printf("  chemcalc list\n");
    printf("  chemcalc show ID\n");
    printf("  chemcalc yield ID ACTUAL_VALUE [--unit g|mg|kg|ug] [--product N]\n");
    printf("  chemcalc rename ID NEW_NAME\n");
    printf("  chemcalc delete ID [-f|--force]\n");
    printf("\nAn ID is either the name you gave the reaction or its random\n");
    printf("hex id; a unique prefix of either is enough.\n");
}

int main(int argc, char **argv) {
    init_colors();
    if (argc < 2) { print_usage(); return 1; }

    if (strcmp(argv[1], "mw") == 0) cmd_mw(argc, argv);
    else if (strcmp(argv[1], "react") == 0) cmd_react();
    else if (strcmp(argv[1], "list") == 0) cmd_list();
    else if (strcmp(argv[1], "show") == 0) cmd_show(argc, argv);
    else if (strcmp(argv[1], "yield") == 0) cmd_yield(argc, argv);
    else if (strcmp(argv[1], "rename") == 0) cmd_rename(argc, argv);
    else if (strcmp(argv[1], "delete") == 0 || strcmp(argv[1], "rm") == 0)
        cmd_delete(argc, argv);
    else if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) print_usage();
    else { fprintf(stderr, "Unknown command '%s'\n\n", argv[1]); print_usage(); return 1; }

    return 0;
}

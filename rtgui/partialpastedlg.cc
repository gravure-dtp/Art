/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "partialpastedlg.h"
#include "multilangmgr.h"
#include "paramsedited.h"
#include "guiutils.h"

PartialPasteDlg::PartialPasteDlg(const Glib::ustring &title, Gtk::Window *parent):
    Gtk::Dialog(title, *parent, true)
{
    set_default_size(700, 600);

    everything_ = Gtk::manage(new Gtk::CheckButton(M("PARTIALPASTE_EVERYTHING")));
    everything_->set_name("PartialPasteHeader");
    everything_conn_ = everything_->signal_toggled().connect(sigc::bind(sigc::mem_fun(*this, &PartialPasteDlg::toggled), everything_));

    Gtk::VBox *vboxes[9];
    Gtk::HSeparator *hseps[9];

    for (int i = 0; i < 9; i++) {
        vboxes[i] = Gtk::manage(new Gtk::VBox());
        vboxes[i]->set_name("PartialPasteGroupContainer");
        hseps[i] = Gtk::manage(new Gtk::HSeparator());
        hseps[i]->set_name("PartialPasteHeaderSep");
    }

    Gtk::VBox *vbCol1 = Gtk::manage(new Gtk::VBox());
    Gtk::VBox *vbCol2 = Gtk::manage(new Gtk::VBox());
    Gtk::VBox *vbCol3 = Gtk::manage(new Gtk::VBox());

    for (int i = 0; i < 3; i++) {
        vbCol1->pack_start(*vboxes[i], Gtk::PACK_SHRINK, 2);
    }

    for (int i = 3; i < 8; i++) {
        vbCol2->pack_start(*vboxes[i], Gtk::PACK_SHRINK, 2);
    }

    for (int i = 8; i < 9; i++) {
        vbCol3->pack_start(*vboxes[i], Gtk::PACK_SHRINK, 2);
    }

    Gtk::VBox *vbtop = Gtk::manage(new Gtk::VBox());
    vbtop->pack_start(*everything_, Gtk::PACK_SHRINK, 2);

    Gtk::Dialog::get_content_area()->pack_start(*vbtop, Gtk::PACK_SHRINK, 2);

    Gtk::HBox *hbmain = Gtk::manage(new Gtk::HBox());
    hbmain->pack_start(*vbCol1);
    Gtk::VSeparator *vsep1 = Gtk::manage(new Gtk::VSeparator());
    setExpandAlignProperties(vsep1, false, true, Gtk::ALIGN_CENTER, Gtk::ALIGN_FILL);
    hbmain->pack_start(*vsep1);
    hbmain->pack_start(*vbCol2);
    Gtk::VSeparator *vsep2 = Gtk::manage(new Gtk::VSeparator());
    setExpandAlignProperties(vsep2, false, true, Gtk::ALIGN_CENTER, Gtk::ALIGN_FILL);
    hbmain->pack_start(*vsep2);
    hbmain->pack_start(*vbCol3);

    Gtk::ScrolledWindow *scrolledwindow = Gtk::manage(new Gtk::ScrolledWindow());
    scrolledwindow->set_name("PartialPaste");
    scrolledwindow->set_can_focus(true);
    scrolledwindow->set_shadow_type(Gtk::SHADOW_NONE);
    scrolledwindow->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);
    scrolledwindow->property_window_placement().set_value(Gtk::CORNER_TOP_LEFT);

    scrolledwindow->add(*hbmain);

    Gtk::Dialog::get_content_area()->pack_start(*scrolledwindow, Gtk::PACK_EXPAND_WIDGET, 2);

    static const std::vector<std::string> groups = {
        "PARTIALPASTE_BASICGROUP",
        "PARTIALPASTE_DETAILGROUP",
        "PARTIALPASTE_COLORGROUP",
        "PARTIALPASTE_LENSGROUP",
        "PARTIALPASTE_COMPOSITIONGROUP",
        "PARTIALPASTE_LOCALGROUP",
        "PARTIALPASTE_METAGROUP",
        "PARTIALPASTE_RAWGROUP"
    };
    std::vector<Gtk::CheckButton *> gbtns;
    size_t i = 0;
    for (const auto &g : groups) {
        auto b = Gtk::manage(new Gtk::CheckButton(M(g)));
        gbtns.push_back(b);
        b->set_name("PartialPasteHeader");
        vboxes[i]->pack_start(*b, Gtk::PACK_SHRINK, 2);
        vboxes[i]->pack_start(*hseps[i], Gtk::PACK_SHRINK, 2);
        ++i;

        auto conn = b->signal_toggled().connect(sigc::bind(sigc::mem_fun(*this, &PartialPasteDlg::toggled), b));
        buttons_[b] = { conn, {}, true, nullptr };
    }

    struct ToggleInfo {
        std::string label;
        bool *edited;
        int index;
    };
    const std::vector<ToggleInfo> toggles = {
        {"PARTIALPASTE_WHITEBALANCE", &pedited_.wb, 0},
        {"PARTIALPASTE_EXPOSURE", &pedited_.toneCurve, 0},
        {"PARTIALPASTE_TM_LOG", &pedited_.logenc, 0},
        {"PARTIALPASTE_SHADOWSHIGHLIGHTS", &pedited_.sh, 0},
        {"PARTIALPASTE_TONE_EQUALIZER", &pedited_.toneEqualizer, 0},
        {"PARTIALPASTE_TM_FATTAL", &pedited_.fattal, 0},
        {"PARTIALPASTE_PCVIGNETTE", &pedited_.pcvignette, 0},
        {"PARTIALPASTE_GRADIENT", &pedited_.gradient, 0},
        {"PARTIALPASTE_CHANNELMIXERBW", &pedited_.blackwhite, 0},

        {"PARTIALPASTE_SHARPENING", &pedited_.sharpening, 1},
        {"PARTIALPASTE_LOCALCONTRAST", &pedited_.localContrast, 1},
        {"PARTIALPASTE_SHARPENMICRO", &pedited_.sharpenMicro, 1},
        {"PARTIALPASTE_IMPULSEDENOISE", &pedited_.impulseDenoise, 1},
        {"PARTIALPASTE_DIRPYRDENOISE", &pedited_.denoise, 1},
        {"PARTIALPASTE_DEFRINGE", &pedited_.defringe, 1},
        {"PARTIALPASTE_DEHAZE", &pedited_.dehaze, 1},
        {"PARTIALPASTE_GRAIN", &pedited_.grain, 1},

        {"PARTIALPASTE_ICMSETTINGS", &pedited_.icm, 2},
        {"PARTIALPASTE_CHANNELMIXER", &pedited_.chmixer, 2},
        {"PARTIALPASTE_FILMSIMULATION", &pedited_.filmSimulation, 2},
        {"PARTIALPASTE_SOFTLIGHT", &pedited_.softlight, 2},
        {"PARTIALPASTE_RGBCURVES", &pedited_.rgbCurves, 2},
        {"PARTIALPASTE_LABCURVE", &pedited_.labCurve, 2},

        {"PARTIALPASTE_DISTORTION", &pedited_.distortion, 3},
        {"PARTIALPASTE_CACORRECTION", &pedited_.cacorrection, 3},
        {"PARTIALPASTE_VIGNETTING", &pedited_.vignetting, 3},
        {"PARTIALPASTE_LENSPROFILE", &pedited_.lensProf, 3},

        {"PARTIALPASTE_COARSETRANS", &pedited_.coarse, 4},
        {"PARTIALPASTE_ROTATION", &pedited_.rotate, 4},
        {"PARTIALPASTE_CROP", &pedited_.crop, 4},
        {"PARTIALPASTE_RESIZE", &pedited_.resize, 4},
        {"PARTIALPASTE_PRSHARPENING", &pedited_.prsharpening, 4},
        {"PARTIALPASTE_PERSPECTIVE", &pedited_.perspective, 4},
        {"PARTIALPASTE_COMMONTRANSFORMPARAMS", &pedited_.commonTrans, 4},
        
        {"PARTIALPASTE_COLORCORRECTION", &pedited_.colorcorrection, 5},
        {"PARTIALPASTE_DIRPYREQUALIZER", &pedited_.dirpyrequalizer, 5},
        {"PARTIALPASTE_SMOOTHING", &pedited_.smoothing, 5},
        {"PARTIALPASTE_EPD", &pedited_.epd, 5},

        {"PARTIALPASTE_METADATA", &pedited_.metadata, 6},
        {"PARTIALPASTE_EXIFCHANGES", &pedited_.exif, 6},
        {"PARTIALPASTE_IPTCINFO", &pedited_.iptc, 6},

        {"PARTIALPASTE_RAW_DEMOSAIC", &pedited_.demosaic, 7},
        {"PARTIALPASTE_RAW_BLACK", &pedited_.rawBlack, 7},
        {"PARTIALPASTE_RAW_WHITE", &pedited_.rawWhite, 7},
        {"PARTIALPASTE_RAW_PREPROCESSING", &pedited_.rawPreprocessing, 7},
        {"PARTIALPASTE_HOT_DEAD_PIXEL_FILTER", &pedited_.hotDeadPixelFilter, 7},
        {"PARTIALPASTE_DARKFRAME", &pedited_.darkframe, 7},
        {"PARTIALPASTE_FLATFIELD", &pedited_.flatfield, 7},
        {"PARTIALPASTE_RAWCA", &pedited_.rawCA, 7}
    };

    for (const auto &t : toggles) {
        auto b = Gtk::manage(new Gtk::CheckButton(M(t.label)));
        auto conn = b->signal_toggled().connect(sigc::bind(sigc::mem_fun(*this, &PartialPasteDlg::toggled), b));
        int i = t.index;
        auto master = gbtns[i];
        buttons_[b] = {conn, {master}, false, t.edited};
        buttons_[master].related.push_back(b);
        vboxes[i]->pack_start(*b, Gtk::PACK_SHRINK, 2);
    }
    
    hbmain->show();
    scrolledwindow->show();
    

    // basic       = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_BASICGROUP")));
    // basic       ->set_name("PartialPasteHeader");
    // detail      = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DETAILGROUP")));
    // detail      ->set_name("PartialPasteHeader");
    // color       = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_COLORGROUP")));
    // color       ->set_name("PartialPasteHeader");
    // lens        = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_LENSGROUP")));
    // lens        ->set_name("PartialPasteHeader");
    // composition = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_COMPOSITIONGROUP")));
    // composition ->set_name("PartialPasteHeader");
    // meta        = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_METAGROUP")));
    // meta        ->set_name("PartialPasteHeader");
    // raw         = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAWGROUP")));
    // raw         ->set_name("PartialPasteHeader");
    // local = Gtk::manage(new Gtk::CheckButton(M("PARTIALPASTE_LOCALGROUP")));
    // local->set_name("PartialPasteHeader");
    

    // // Basic Settings:
    // wb          = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_WHITEBALANCE")));
    // exposure    = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_EXPOSURE")));
    // logenc      = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_TM_LOG")));
    // sh          = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_SHADOWSHIGHLIGHTS")));
    // toneEqualizer = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_TONE_EQUALIZER")));
    // fattal      = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_TM_FATTAL")));
    // pcvignette  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_PCVIGNETTE")));
    // gradient    = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_GRADIENT")));
    // labcurve    = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_LABCURVE")));

    // // Detail Settings:
    // sharpen     = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_SHARPENING")));
    // localcontrast = Gtk::manage(new Gtk::CheckButton(M("PARTIALPASTE_LOCALCONTRAST")));
    // sharpenmicro = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_SHARPENMICRO")));
    // impden      = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_IMPULSEDENOISE")));
    // dirpyrden   = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DIRPYRDENOISE")));
    // defringe    = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DEFRINGE")));
    // //dirpyreq    = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DIRPYREQUALIZER")));
    // dehaze = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DEHAZE")) );
    // grain = Gtk::manage(new Gtk::CheckButton(M("PARTIALPASTE_GRAIN")));

    // // Color-Related Settings
    // icm         = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_ICMSETTINGS")));
    // chmixer     = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_CHANNELMIXER")));
    // blackwhite  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_CHANNELMIXERBW")));
    // filmSimulation = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_FILMSIMULATION")) );
    // softlight = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_SOFTLIGHT")) );
    // rgbcurves   = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RGBCURVES")));

    // // Lens-Related Settings
    // distortion  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DISTORTION")));
    // cacorr      = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_CACORRECTION")));
    // vignetting  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_VIGNETTING")));
    // lcp         = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_LENSPROFILE")));

    // // Composition Settings:
    // coarserot    = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_COARSETRANS")));
    // finerot      = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_ROTATION")));
    // crop         = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_CROP")));
    // resize       = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RESIZE")));
    // prsharpening = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_PRSHARPENING")));
    // perspective  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_PERSPECTIVE")));
    // commonTrans  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_COMMONTRANSFORMPARAMS")));

    // // Metadata:
    // metadata = Gtk::manage(new Gtk::CheckButton(M("PARTIALPASTE_METADATA")));
    // exifch      = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_EXIFCHANGES")));
    // iptc        = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_IPTCINFO")));

    // // Local settings:
    // colorcorrection = Gtk::manage(new Gtk::CheckButton(M("PARTIALPASTE_COLORCORRECTION")));
    // cbdl = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DIRPYREQUALIZER")));
    // smoothing = Gtk::manage(new Gtk::CheckButton(M("PARTIALPASTE_SMOOTHING")));
    // epd         = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_EPD")));

    // // Raw Settings:
    // raw_method          = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAW_DMETHOD")));
    // raw_imagenum        = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAW_IMAGENUM")));
    // raw_border          = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAW_BORDER")));
    // raw_pixelshift      = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAW_PIXELSHIFT")));
    // raw_ccSteps         = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAW_FALSECOLOR")));
    // raw_dcb_iterations  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAW_DCBITERATIONS")));
    // raw_dcb_enhance     = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAW_DCBENHANCE")));
    // raw_lmmse_iterations = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAW_LMMSEITERATIONS")));
    // //---
    // raw_linenoise       = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_PREPROCESS_LINEDENOISE")));
    // raw_greenthresh     = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_PREPROCESS_GREENEQUIL")));
    // raw_hotpix_filt     = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_PREPROCESS_HOTPIXFILT")));
    // raw_deadpix_filt    = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_PREPROCESS_DEADPIXFILT")));
    // raw_pdaf_lines_filter = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_PREPROCESS_PDAFLINESFILTER")));
    // //---
    // raw_expos           = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAWEXPOS_LINEAR")));
    // raw_black           = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAWEXPOS_BLACK")));
    // //---
    // df_file             = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DARKFRAMEFILE")));
    // df_AutoSelect       = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DARKFRAMEAUTOSELECT")));
    // //---
    // ff_file             = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_FLATFIELDFILE")));
    // ff_AutoSelect       = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_FLATFIELDAUTOSELECT")));
    // ff_BlurType         = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_FLATFIELDBLURTYPE")));
    // ff_BlurRadius       = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_FLATFIELDBLURRADIUS")));
    // ff_ClipControl      = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_FLATFIELDCLIPCONTROL")));
    // //---
    // raw_ca_autocorrect  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAWCACORR_AUTO")));
    // raw_caredblue       = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAWCACORR_CAREDBLUE")));
    // raw_ca_avoid_colourshift = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAWCACORR_AVOIDCOLORSHIFT")));

    // Gtk::VBox* vboxes[9];
    // Gtk::HSeparator* hseps[9];

    // for (int i = 0; i < 9; i++) {
    //     vboxes[i] = Gtk::manage (new Gtk::VBox ());
    //     vboxes[i]->set_name("PartialPasteGroupContainer");
    //     hseps[i] = Gtk::manage (new Gtk::HSeparator ());
    //     hseps[i]->set_name("PartialPasteHeaderSep");
    // }

    // //BASIC
    // vboxes[0]->pack_start (*basic, Gtk::PACK_SHRINK, 2);
    // vboxes[0]->pack_start (*hseps[0], Gtk::PACK_SHRINK, 2);
    // vboxes[0]->pack_start (*wb, Gtk::PACK_SHRINK, 2);
    // vboxes[0]->pack_start (*exposure, Gtk::PACK_SHRINK, 2);
    // vboxes[0]->pack_start (*logenc, Gtk::PACK_SHRINK, 2);
    // vboxes[0]->pack_start (*sh, Gtk::PACK_SHRINK, 2);
    // //vboxes[0]->pack_start (*epd, Gtk::PACK_SHRINK, 2);
    // vboxes[0]->pack_start (*fattal, Gtk::PACK_SHRINK, 2);
    // vboxes[0]->pack_start (*toneEqualizer, Gtk::PACK_SHRINK, 2);
    // vboxes[0]->pack_start (*pcvignette, Gtk::PACK_SHRINK, 2);
    // vboxes[0]->pack_start (*gradient, Gtk::PACK_SHRINK, 2);
    // //vboxes[0]->pack_start (*labcurve, Gtk::PACK_SHRINK, 2);
    // vboxes[0]->pack_start(*blackwhite, Gtk::PACK_SHRINK, 2);

    // //DETAIL
    // vboxes[1]->pack_start (*detail, Gtk::PACK_SHRINK, 2);
    // vboxes[1]->pack_start (*hseps[1], Gtk::PACK_SHRINK, 2);
    // vboxes[1]->pack_start (*sharpen, Gtk::PACK_SHRINK, 2);
    // vboxes[1]->pack_start (*localcontrast, Gtk::PACK_SHRINK, 2);    
    // vboxes[1]->pack_start (*sharpenmicro, Gtk::PACK_SHRINK, 2);
    // vboxes[1]->pack_start (*impden, Gtk::PACK_SHRINK, 2);
    // vboxes[1]->pack_start (*dirpyrden, Gtk::PACK_SHRINK, 2);
    // vboxes[1]->pack_start (*defringe, Gtk::PACK_SHRINK, 2);
    // //vboxes[1]->pack_start (*dirpyreq, Gtk::PACK_SHRINK, 2);
    // vboxes[1]->pack_start (*dehaze, Gtk::PACK_SHRINK, 2);
    // vboxes[1]->pack_start (*grain, Gtk::PACK_SHRINK, 2);
    // //vboxes[1]->pack_start(*smoothing, Gtk::PACK_SHRINK, 2);

    // //COLOR
    // vboxes[2]->pack_start (*color, Gtk::PACK_SHRINK, 2);
    // vboxes[2]->pack_start (*hseps[2], Gtk::PACK_SHRINK, 2);
    // vboxes[2]->pack_start (*icm, Gtk::PACK_SHRINK, 2);
    // vboxes[2]->pack_start (*chmixer, Gtk::PACK_SHRINK, 2);
    // //vboxes[2]->pack_start (*blackwhite, Gtk::PACK_SHRINK, 2);
    // vboxes[2]->pack_start (*filmSimulation, Gtk::PACK_SHRINK, 2);
    // vboxes[2]->pack_start (*softlight, Gtk::PACK_SHRINK, 2);
    // vboxes[2]->pack_start (*rgbcurves, Gtk::PACK_SHRINK, 2);
    // vboxes[2]->pack_start (*labcurve, Gtk::PACK_SHRINK, 2);    

    // //LENS
    // vboxes[3]->pack_start (*lens, Gtk::PACK_SHRINK, 2);
    // vboxes[3]->pack_start (*hseps[3], Gtk::PACK_SHRINK, 2);
    // vboxes[3]->pack_start (*distortion, Gtk::PACK_SHRINK, 2);
    // vboxes[3]->pack_start (*cacorr, Gtk::PACK_SHRINK, 2);
    // vboxes[3]->pack_start (*vignetting, Gtk::PACK_SHRINK, 2);
    // vboxes[3]->pack_start (*lcp, Gtk::PACK_SHRINK, 2);

    // //COMPOSITION
    // vboxes[4]->pack_start (*composition, Gtk::PACK_SHRINK, 2);
    // vboxes[4]->pack_start (*hseps[4], Gtk::PACK_SHRINK, 2);
    // vboxes[4]->pack_start (*coarserot, Gtk::PACK_SHRINK, 2);
    // vboxes[4]->pack_start (*finerot, Gtk::PACK_SHRINK, 2);
    // vboxes[4]->pack_start (*crop, Gtk::PACK_SHRINK, 2);
    // vboxes[4]->pack_start (*resize, Gtk::PACK_SHRINK, 2);
    // vboxes[4]->pack_start (*prsharpening, Gtk::PACK_SHRINK, 2);
    // vboxes[4]->pack_start (*perspective, Gtk::PACK_SHRINK, 2);
    // vboxes[4]->pack_start (*commonTrans, Gtk::PACK_SHRINK, 2);

    // //LOCAL
    // vboxes[5]->pack_start(*local, Gtk::PACK_SHRINK, 2);
    // vboxes[5]->pack_start(*hseps[5], Gtk::PACK_SHRINK, 2);
    // vboxes[5]->pack_start(*colorcorrection, Gtk::PACK_SHRINK, 2);
    // vboxes[5]->pack_start(*smoothing, Gtk::PACK_SHRINK, 2);
    // vboxes[5]->pack_start(*cbdl, Gtk::PACK_SHRINK, 2);
    // vboxes[5]->pack_start(*epd, Gtk::PACK_SHRINK, 2);

    // //META
    // vboxes[7]->pack_start (*meta, Gtk::PACK_SHRINK, 2);
    // vboxes[7]->pack_start (*hseps[7], Gtk::PACK_SHRINK, 2);
    // vboxes[7]->pack_start(*metadata, Gtk::PACK_SHRINK, 2);
    // vboxes[7]->pack_start (*exifch, Gtk::PACK_SHRINK, 2);
    // vboxes[7]->pack_start (*iptc, Gtk::PACK_SHRINK, 2);

    // //RAW
    // vboxes[8]->pack_start (*raw, Gtk::PACK_SHRINK, 2);
    // vboxes[8]->pack_start (*hseps[8], Gtk::PACK_SHRINK, 2);
    // vboxes[8]->pack_start (*raw_method, Gtk::PACK_SHRINK, 2);
    // vboxes[8]->pack_start (*raw_border, Gtk::PACK_SHRINK, 2);
    // vboxes[8]->pack_start (*raw_imagenum, Gtk::PACK_SHRINK, 2);
    // vboxes[8]->pack_start (*raw_pixelshift, Gtk::PACK_SHRINK, 2);
    // vboxes[8]->pack_start (*raw_ccSteps, Gtk::PACK_SHRINK, 2);
    // vboxes[8]->pack_start (*raw_dcb_iterations, Gtk::PACK_SHRINK, 2);
    // vboxes[8]->pack_start (*raw_dcb_enhance, Gtk::PACK_SHRINK, 2);
    // vboxes[8]->pack_start (*raw_lmmse_iterations, Gtk::PACK_SHRINK, 2);
    // vboxes[8]->pack_start (*Gtk::manage (new Gtk::HSeparator ()), Gtk::PACK_SHRINK, 0);
    // vboxes[8]->pack_start (*raw_linenoise, Gtk::PACK_SHRINK, 2);
    // vboxes[8]->pack_start (*raw_greenthresh, Gtk::PACK_SHRINK, 2);
    // vboxes[8]->pack_start (*raw_hotpix_filt, Gtk::PACK_SHRINK, 2);
    // vboxes[8]->pack_start (*raw_deadpix_filt, Gtk::PACK_SHRINK, 2);
    // vboxes[8]->pack_start (*raw_pdaf_lines_filter, Gtk::PACK_SHRINK, 2);
    // vboxes[8]->pack_start (*Gtk::manage (new Gtk::HSeparator ()), Gtk::PACK_SHRINK, 0);
    // vboxes[8]->pack_start (*raw_expos, Gtk::PACK_SHRINK, 2);
    // vboxes[8]->pack_start (*raw_black, Gtk::PACK_SHRINK, 2);
    // vboxes[8]->pack_start (*Gtk::manage (new Gtk::HSeparator ()), Gtk::PACK_SHRINK, 0);
    // vboxes[8]->pack_start (*df_file, Gtk::PACK_SHRINK, 2);
    // vboxes[8]->pack_start (*df_AutoSelect, Gtk::PACK_SHRINK, 2);
    // vboxes[8]->pack_start (*Gtk::manage (new Gtk::HSeparator ()), Gtk::PACK_SHRINK, 0);
    // vboxes[8]->pack_start (*ff_file, Gtk::PACK_SHRINK, 2);
    // vboxes[8]->pack_start (*ff_AutoSelect, Gtk::PACK_SHRINK, 2);
    // vboxes[8]->pack_start (*ff_BlurType, Gtk::PACK_SHRINK, 2);
    // vboxes[8]->pack_start (*ff_BlurRadius, Gtk::PACK_SHRINK, 2);
    // vboxes[8]->pack_start (*ff_ClipControl, Gtk::PACK_SHRINK, 2);
    // vboxes[8]->pack_start (*Gtk::manage (new Gtk::HSeparator ()), Gtk::PACK_SHRINK, 0);
    // vboxes[8]->pack_start (*raw_ca_autocorrect, Gtk::PACK_SHRINK, 2);
    // vboxes[8]->pack_start (*raw_caredblue, Gtk::PACK_SHRINK, 2);
    // vboxes[8]->pack_start (*raw_ca_avoid_colourshift, Gtk::PACK_SHRINK, 2);

    // Gtk::VBox* vbCol1 = Gtk::manage (new Gtk::VBox ());
    // Gtk::VBox* vbCol2 = Gtk::manage (new Gtk::VBox ());
    // Gtk::VBox* vbCol3 = Gtk::manage (new Gtk::VBox ());

    // for (int i = 0; i < 3; i++) {
    //     vbCol1->pack_start (*vboxes[i], Gtk::PACK_SHRINK, 2);
    // }

    // for (int i = 3; i < 8; i++) {
    //     vbCol2->pack_start (*vboxes[i], Gtk::PACK_SHRINK, 2);
    // }

    // for (int i = 8; i < 9; i++) {
    //     vbCol3->pack_start (*vboxes[i], Gtk::PACK_SHRINK, 2);
    // }

    // Gtk::VBox* vbtop = Gtk::manage (new Gtk::VBox ());
    // vbtop->pack_start (*everything, Gtk::PACK_SHRINK, 2);

    // Gtk::Dialog::get_content_area()->pack_start (*vbtop, Gtk::PACK_SHRINK, 2);

    // Gtk::HBox* hbmain = Gtk::manage (new Gtk::HBox ());
    // hbmain->pack_start (*vbCol1);
    // Gtk::VSeparator *vsep1 = Gtk::manage (new Gtk::VSeparator ());
    // setExpandAlignProperties(vsep1, false, true, Gtk::ALIGN_CENTER, Gtk::ALIGN_FILL);
    // hbmain->pack_start (*vsep1);
    // hbmain->pack_start (*vbCol2);
    // Gtk::VSeparator *vsep2 = Gtk::manage (new Gtk::VSeparator ());
    // setExpandAlignProperties(vsep2, false, true, Gtk::ALIGN_CENTER, Gtk::ALIGN_FILL);
    // hbmain->pack_start (*vsep2);
    // hbmain->pack_start (*vbCol3);

    // scrolledwindow = Gtk::manage ( new Gtk::ScrolledWindow() );
    // scrolledwindow->set_name("PartialPaste");
    // scrolledwindow->set_can_focus(true);
    // scrolledwindow->set_shadow_type(Gtk::SHADOW_NONE);
    // scrolledwindow->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);
    // scrolledwindow->property_window_placement().set_value(Gtk::CORNER_TOP_LEFT);

    // scrolledwindow->add(*hbmain);

    // Gtk::Dialog::get_content_area()->pack_start (*scrolledwindow, Gtk::PACK_EXPAND_WIDGET, 2);

    // hbmain->show();
    // scrolledwindow->show ();

    // // This can be improved
    // // there is currently no binding of subsettings to CheckButton 'everything' for its inconsistent status
    // everythingConn  = everything->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::everythingToggled));
    // basicConn       = basic->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::basicToggled));
    // detailConn      = detail->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::detailToggled));
    // colorConn       = color->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::colorToggled));
    // lensConn        = lens->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::lensToggled));
    // compositionConn = composition->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::compositionToggled));
    // metaConn        = meta->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::metaToggled));
    // rawConn         = raw->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::rawToggled));
    // localConn = local->signal_toggled().connect(sigc::mem_fun(*this, &PartialPasteDlg::localToggled));


    // // Basic Settings
    // wbConn          = wb->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));
    // exposureConn    = exposure->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));
    // logencConn      = logenc->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));
    // shConn          = sh->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));
    // toneEqualizerConn = toneEqualizer->signal_toggled().connect(sigc::bind(sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));
    // fattalConn      = fattal->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));
    // pcvignetteConn  = pcvignette->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));
    // gradientConn    = gradient->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));
    // labcurveConn    = labcurve->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));

    // // Detail Settings:
    // sharpenConn     = sharpen->signal_toggled().connect (sigc::bind (sigc::mem_fun(*detail, &Gtk::CheckButton::set_inconsistent), true));
    // localcontrastConn = localcontrast->signal_toggled().connect (sigc::bind (sigc::mem_fun(*detail, &Gtk::CheckButton::set_inconsistent), true));
    // microcontrastConn = sharpenmicro->signal_toggled().connect (sigc::bind (sigc::mem_fun(*detail, &Gtk::CheckButton::set_inconsistent), true));
    // impdenConn      = impden->signal_toggled().connect (sigc::bind (sigc::mem_fun(*detail, &Gtk::CheckButton::set_inconsistent), true));
    // dirpyrdenConn   = dirpyrden->signal_toggled().connect (sigc::bind (sigc::mem_fun(*detail, &Gtk::CheckButton::set_inconsistent), true));
    // defringeConn    = defringe->signal_toggled().connect (sigc::bind (sigc::mem_fun(*detail, &Gtk::CheckButton::set_inconsistent), true));
    // //dirpyreqConn    = dirpyreq->signal_toggled().connect (sigc::bind (sigc::mem_fun(*detail, &Gtk::CheckButton::set_inconsistent), true));
    // dehazeConn    = dehaze->signal_toggled().connect (sigc::bind (sigc::mem_fun(*detail, &Gtk::CheckButton::set_inconsistent), true));
    // grainConn    = grain->signal_toggled().connect (sigc::bind (sigc::mem_fun(*detail, &Gtk::CheckButton::set_inconsistent), true));

    // // Color-related Settings:
    // icmConn         = icm->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));
    // chmixerConn     = chmixer->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));
    // chmixerbwConn   = blackwhite->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));
    // filmSimulationConn = filmSimulation->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));
    // softlightConn = softlight->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));
    // rgbcurvesConn   = rgbcurves->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));

    // // Lens-Related Settings:
    // distortionConn  = distortion->signal_toggled().connect (sigc::bind (sigc::mem_fun(*lens, &Gtk::CheckButton::set_inconsistent), true));
    // cacorrConn      = cacorr->signal_toggled().connect (sigc::bind (sigc::mem_fun(*lens, &Gtk::CheckButton::set_inconsistent), true));
    // vignettingConn  = vignetting->signal_toggled().connect (sigc::bind (sigc::mem_fun(*lens, &Gtk::CheckButton::set_inconsistent), true));
    // lcpConn         = lcp->signal_toggled().connect (sigc::bind (sigc::mem_fun(*lens, &Gtk::CheckButton::set_inconsistent), true));

    // // Composition Settings:
    // coarserotConn   = coarserot->signal_toggled().connect (sigc::bind (sigc::mem_fun(*composition, &Gtk::CheckButton::set_inconsistent), true));
    // finerotConn     = finerot->signal_toggled().connect (sigc::bind (sigc::mem_fun(*composition, &Gtk::CheckButton::set_inconsistent), true));
    // cropConn        = crop->signal_toggled().connect (sigc::bind (sigc::mem_fun(*composition, &Gtk::CheckButton::set_inconsistent), true));
    // resizeConn      = resize->signal_toggled().connect (sigc::bind (sigc::mem_fun(*composition, &Gtk::CheckButton::set_inconsistent), true));
    // prsharpeningConn = prsharpening->signal_toggled().connect (sigc::bind (sigc::mem_fun(*composition, &Gtk::CheckButton::set_inconsistent), true));
    // perspectiveConn = perspective->signal_toggled().connect (sigc::bind (sigc::mem_fun(*composition, &Gtk::CheckButton::set_inconsistent), true));
    // commonTransConn = commonTrans->signal_toggled().connect (sigc::bind (sigc::mem_fun(*composition, &Gtk::CheckButton::set_inconsistent), true));

    // // Metadata:
    // metadataConn = metadata->signal_toggled().connect(sigc::bind (sigc::mem_fun(*meta, &Gtk::CheckButton::set_inconsistent), true));
    // exifchConn      = exifch->signal_toggled().connect (sigc::bind (sigc::mem_fun(*meta, &Gtk::CheckButton::set_inconsistent), true));
    // iptcConn        = iptc->signal_toggled().connect (sigc::bind (sigc::mem_fun(*meta, &Gtk::CheckButton::set_inconsistent), true));

    // // Local:
    // colorcorrectionConn = colorcorrection->signal_toggled().connect(sigc::bind(sigc::mem_fun(*local, &Gtk::CheckButton::set_inconsistent), true));
    // smoothingConn = smoothing->signal_toggled().connect(sigc::bind(sigc::mem_fun(*local, &Gtk::CheckButton::set_inconsistent), true));
    // cbdlConn = cbdl->signal_toggled().connect(sigc::bind(sigc::mem_fun(*local, &Gtk::CheckButton::set_inconsistent), true));
    // epdConn         = epd->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));

    // // Raw Settings:
    // raw_methodConn          = raw_method->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    // raw_borderConn          = raw_border->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    // raw_imagenumConn        = raw_imagenum->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    // raw_pixelshiftConn      = raw_pixelshift->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    // raw_ccStepsConn         = raw_ccSteps->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    // raw_dcb_iterationsConn  = raw_dcb_iterations->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    // raw_dcb_enhanceConn     = raw_dcb_enhance->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    // raw_lmmse_iterationsConn  = raw_lmmse_iterations->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    // //---
    // raw_linenoiseConn       = raw_linenoise->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    // raw_greenthreshConn     = raw_greenthresh->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    // raw_hotpix_filtConn     = raw_hotpix_filt->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    // raw_deadpix_filtConn    = raw_deadpix_filt->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    // raw_pdaf_lines_filterConn = raw_pdaf_lines_filter->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    // //---
    // raw_exposConn           = raw_expos->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    // raw_blackConn           = raw_black->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    // //---
    // df_fileConn             = df_file->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    // df_AutoSelectConn       = df_AutoSelect->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    // //---
    // ff_fileConn             = ff_file->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    // ff_AutoSelectConn       = ff_AutoSelect->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    // ff_BlurTypeConn         = ff_BlurType->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    // ff_BlurRadiusConn       = ff_BlurRadius->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    // ff_ClipControlConn      = ff_ClipControl->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    // //---
    // raw_ca_autocorrectConn  = raw_ca_autocorrect->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    // raw_caredblueConn       = raw_caredblue->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    // raw_ca_avoid_colourshiftconn = raw_ca_avoid_colourshift->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));

    add_button(M("GENERAL_OK"), Gtk::RESPONSE_OK);
    add_button(M("GENERAL_CANCEL"), Gtk::RESPONSE_CANCEL);
    set_response_sensitive(Gtk::RESPONSE_OK);
    set_default_response(Gtk::RESPONSE_OK);
    show_all_children();
}


void PartialPasteDlg::toggled(Gtk::CheckButton *which)
{
    which->set_inconsistent(false);
    
    if (which == everything_) {
        for (auto &p : buttons_) {
            auto b = p.first;
            ConnectionBlocker blocker(p.second.conn);
            b->set_active(which->get_active());
        }
    } else {
        auto &g = buttons_[which];
        if (!which->get_active()) {
            ConnectionBlocker blocker(everything_conn_);
            everything_->set_active(false);
        }
        if (g.is_master || !which->get_active()) {
            for (auto c : g.related) {
                auto &ci = buttons_[c];
                ConnectionBlocker blocker(ci.conn);
                c->set_active(which->get_active());
            }
        }
    }
}


ParamsEdited PartialPasteDlg::getParamsEdited()
{
    for (auto &p : buttons_) {
        if (p.second.edited) {
            *p.second.edited = p.first->get_active();
        }
    }
    return pedited_;
}

// void PartialPasteDlg::everythingToggled ()
// {

//     ConnectionBlocker basicBlocker(basicConn);
//     ConnectionBlocker detailBlocker(detailConn);
//     ConnectionBlocker colorBlocker(colorConn);
//     ConnectionBlocker lensBlocker(lensConn);
//     ConnectionBlocker compositionBlocker(compositionConn);
//     ConnectionBlocker metaBlocker(metaConn);
//     ConnectionBlocker rawBlocker(rawConn);
//     ConnectionBlocker localBlocker(localConn);

//     everything->set_inconsistent (false);

//     //toggle group headings
//     basic->set_active(everything->get_active());
//     detail->set_active(everything->get_active());
//     color->set_active(everything->get_active());
//     lens->set_active(everything->get_active());
//     composition->set_active(everything->get_active());
//     meta->set_active(everything->get_active());
//     raw->set_active(everything->get_active());
//     local->set_active(everything->get_active());

//     //toggle group children
//     PartialPasteDlg::basicToggled ();
//     PartialPasteDlg::detailToggled ();
//     PartialPasteDlg::colorToggled ();
//     PartialPasteDlg::lensToggled ();
//     PartialPasteDlg::compositionToggled ();
//     PartialPasteDlg::metaToggled ();
//     PartialPasteDlg::rawToggled ();
//     PartialPasteDlg::localToggled();
// }

// void PartialPasteDlg::rawToggled ()
// {

//     ConnectionBlocker raw_methodBlocker(raw_methodConn);
//     ConnectionBlocker raw_borderBlocker(raw_borderConn);
//     ConnectionBlocker raw_imagenumBlocker(raw_imagenumConn);
//     ConnectionBlocker raw_pixelshiftBlocker(raw_pixelshiftConn);
//     ConnectionBlocker raw_ccStepsBlocker(raw_ccStepsConn);
//     ConnectionBlocker raw_dcb_iterationsBlocker(raw_dcb_iterationsConn);
//     ConnectionBlocker raw_dcb_enhanceBlocker(raw_dcb_enhanceConn);
//     ConnectionBlocker raw_lmmse_iterationsBlocker(raw_lmmse_iterationsConn);
//     ConnectionBlocker raw_linenoiseBlocker(raw_linenoiseConn);
//     ConnectionBlocker raw_greenthreshBlocker(raw_greenthreshConn);
//     ConnectionBlocker raw_hotpix_filtBlocker(raw_hotpix_filtConn);
//     ConnectionBlocker raw_deadpix_filtBlocker(raw_deadpix_filtConn);
//     ConnectionBlocker raw_pdaf_lines_filterBlocker(raw_pdaf_lines_filterConn);
//     ConnectionBlocker raw_exposBlocker(raw_exposConn);
//     ConnectionBlocker raw_blackBlocker(raw_blackConn);
//     ConnectionBlocker df_fileBlocker(df_fileConn);
//     ConnectionBlocker df_AutoSelectBlocker(df_AutoSelectConn);
//     ConnectionBlocker ff_fileBlocker(ff_fileConn);
//     ConnectionBlocker ff_AutoSelectBlocker(ff_AutoSelectConn);
//     ConnectionBlocker ff_BlurTypeBlocker(ff_BlurTypeConn);
//     ConnectionBlocker ff_BlurRadiusBlocker(ff_BlurRadiusConn);
//     ConnectionBlocker ff_ClipControlBlocker(ff_ClipControlConn);
//     ConnectionBlocker raw_ca_autocorrectBlocker(raw_ca_autocorrectConn);
//     ConnectionBlocker raw_caredblueBlocker(raw_caredblueConn);
//     ConnectionBlocker raw_ca_avoid_colourshiftBlocker(raw_ca_avoid_colourshiftconn);

//     raw->set_inconsistent (false);

//     raw_method->set_active (raw->get_active ());
//     raw_border->set_active (raw->get_active ());
//     raw_imagenum->set_active (raw->get_active ());
//     raw_pixelshift->set_active (raw->get_active ());
//     raw_ccSteps->set_active (raw->get_active ());
//     raw_dcb_iterations->set_active (raw->get_active ());
//     raw_dcb_enhance->set_active (raw->get_active ());
//     raw_lmmse_iterations->set_active (raw->get_active ());
//     raw_linenoise->set_active (raw->get_active ());
//     raw_greenthresh->set_active (raw->get_active ());
//     raw_hotpix_filt->set_active (raw->get_active ());
//     raw_deadpix_filt->set_active (raw->get_active ());
//     raw_pdaf_lines_filter->set_active (raw->get_active ());
//     raw_expos->set_active (raw->get_active ());
//     raw_black->set_active (raw->get_active ());
//     df_file->set_active (raw->get_active ());
//     df_AutoSelect->set_active (raw->get_active ());
//     ff_file->set_active (raw->get_active ());
//     ff_AutoSelect->set_active (raw->get_active ());
//     ff_BlurType->set_active (raw->get_active ());
//     ff_BlurRadius->set_active (raw->get_active ());
//     ff_ClipControl->set_active (raw->get_active ());
//     raw_ca_autocorrect->set_active (raw->get_active ());
//     raw_caredblue->set_active (raw->get_active ());
//     raw_ca_avoid_colourshift->set_active (raw->get_active ());
// }

// void PartialPasteDlg::basicToggled ()
// {

//     ConnectionBlocker wbBlocker(wbConn);
//     ConnectionBlocker exposureBlocker(exposureConn);
//     ConnectionBlocker logencBlocker(logencConn);
//     ConnectionBlocker shBlocker(shConn);
//     ConnectionBlocker toneEqualizerBlocker(toneEqualizerConn);
//     //ConnectionBlocker epdBlocker(epdConn);
//     ConnectionBlocker fattalBlocker(fattalConn);
//     ConnectionBlocker pcvignetteBlocker(pcvignetteConn);
//     ConnectionBlocker gradientBlocker(gradientConn);
//     //ConnectionBlocker labcurveBlocker(labcurveConn);
//     ConnectionBlocker blackwhiteBlocker(chmixerbwConn);

//     basic->set_inconsistent (false);

//     wb->set_active (basic->get_active ());
//     exposure->set_active (basic->get_active ());
//     logenc->set_active (basic->get_active ());
//     sh->set_active (basic->get_active ());
//     toneEqualizer->set_active(basic->get_active());
//     //epd->set_active (basic->get_active ());
//     fattal->set_active (basic->get_active ());
//     pcvignette->set_active (basic->get_active ());
//     gradient->set_active (basic->get_active ());
//     //labcurve->set_active (basic->get_active ());
//     blackwhite->set_active(basic->get_active());
// }

// void PartialPasteDlg::detailToggled ()
// {

//     ConnectionBlocker sharpenBlocker(sharpenConn);
//     ConnectionBlocker localcontrastBlocker(localcontrastConn);
//     ConnectionBlocker microcontrastBlocker(microcontrastConn);
//     ConnectionBlocker impdenBlocker(impdenConn);
//     ConnectionBlocker dirpyrdenBlocker(dirpyrdenConn);
//     ConnectionBlocker defringeBlocker(defringeConn);
//     //ConnectionBlocker dirpyreqBlocker(dirpyreqConn);
//     ConnectionBlocker dehazeBlocker(dehazeConn);
//     ConnectionBlocker grainBlocker(grainConn);
//     //ConnectionBlocker smoothingBlocker(smoothingConn);

//     detail->set_inconsistent (false);

//     sharpen->set_active (detail->get_active ());
//     localcontrast->set_active(detail->get_active());
//     sharpenmicro->set_active (detail->get_active ());
//     impden->set_active (detail->get_active ());
//     dirpyrden->set_active (detail->get_active ());
//     defringe->set_active (detail->get_active ());
//     //dirpyreq->set_active (detail->get_active ());
//     dehaze->set_active (detail->get_active ());
//     grain->set_active (detail->get_active ());
//     //smoothing->set_active(detail->get_active());
// }


// void PartialPasteDlg::colorToggled ()
// {

//     ConnectionBlocker icmBlocker(icmConn);
//     ConnectionBlocker chmixerBlocker(chmixerConn);
//     //ConnectionBlocker chmixerbwBlocker(chmixerbwConn);
//     ConnectionBlocker filmSimulationBlocker(filmSimulationConn);
//     ConnectionBlocker softlightBlocker(softlightConn);
//     ConnectionBlocker rgbcurvesBlocker(rgbcurvesConn);
//     ConnectionBlocker labcurveBlocker(labcurveConn);
    
//     color->set_inconsistent (false);

//     icm->set_active (color->get_active ());
//     chmixer->set_active (color->get_active ());
//     //blackwhite->set_active (color->get_active ());
//     labcurve->set_active(color->get_active());
//     filmSimulation->set_active (color->get_active ());
//     softlight->set_active (color->get_active ());
//     rgbcurves->set_active (color->get_active ());
// }

// void PartialPasteDlg::lensToggled ()
// {

//     ConnectionBlocker distortionBlocker(distortionConn);
//     ConnectionBlocker cacorrBlocker(cacorrConn);
//     ConnectionBlocker vignettingBlocker(vignettingConn);
//     ConnectionBlocker lcpBlocker(lcpConn);

//     lens->set_inconsistent (false);

//     distortion->set_active (lens->get_active ());
//     cacorr->set_active (lens->get_active ());
//     vignetting->set_active (lens->get_active ());
//     lcp->set_active (lens->get_active ());
// }

// void PartialPasteDlg::compositionToggled ()
// {

//     ConnectionBlocker coarserotBlocker(coarserotConn);
//     ConnectionBlocker finerotBlocker(finerotConn);
//     ConnectionBlocker cropBlocker(cropConn);
//     ConnectionBlocker resizeBlocker(resizeConn);
//     ConnectionBlocker prsharpeningBlocker(prsharpeningConn);
//     ConnectionBlocker perspectiveBlocker(perspectiveConn);
//     ConnectionBlocker commonTransBlocker(commonTransConn);

//     composition->set_inconsistent (false);

//     coarserot->set_active (composition->get_active ());
//     finerot->set_active (composition->get_active ());
//     crop->set_active (composition->get_active ());
//     resize->set_active (composition->get_active ());
//     prsharpening->set_active (composition->get_active ());
//     perspective->set_active (composition->get_active ());
//     commonTrans->set_active (composition->get_active ());
// }

// void PartialPasteDlg::metaToggled ()
// {

//     ConnectionBlocker metadataBlocker(metadataConn);
//     ConnectionBlocker exifchBlocker(exifchConn);
//     ConnectionBlocker iptcBlocker(iptcConn);

//     meta->set_inconsistent (false);

//     metadata->set_active(meta->get_active());
//     exifch->set_active (meta->get_active ());
//     iptc->set_active (meta->get_active ());
// }


// void PartialPasteDlg::localToggled()
// {

//     ConnectionBlocker colorcorrectionBlocker(colorcorrectionConn);
//     ConnectionBlocker smoothingBlocker(smoothingConn);
//     ConnectionBlocker cbdlBlocker(cbdlConn);
//     ConnectionBlocker epdBlocker(epdConn);    

//     local->set_inconsistent(false);

//     colorcorrection->set_active(local->get_active());
//     smoothing->set_active(local->get_active());
//     cbdl->set_active(local->get_active());
//     epd->set_active (local->get_active ());
// }


// /*
//  * Copies the selected items from the source ProcParams+ParamsEdited(optional)
//  * to the destination ProcParams.
//  */
// void PartialPasteDlg::applyPaste (rtengine::procparams::ProcParams* dstPP, ParamsEdited* dstPE, const rtengine::procparams::ProcParams* srcPP, const ParamsEdited* srcPE)
// {

//     ParamsEdited falsePE;  // falsePE is a workaround to set a group of ParamsEdited to false
//     ParamsEdited filterPE(true); // Contains the initial information about the loaded values

//     if (srcPE) {
//         filterPE = *srcPE;
//     }

//     // the general section is always ignored, whichever operation we use the PartialPaste for
//     filterPE.general = falsePE.general;

//     // Now we filter out the filter depending on the checked items
//     if (!wb->get_active ()) {
//         filterPE.wb         = falsePE.wb;
//     }

//     if (!exposure->get_active ()) {
//         filterPE.toneCurve  = falsePE.toneCurve;
//     }

//     if (!localcontrast->get_active()) {
//         filterPE.localContrast = falsePE.localContrast;
//     }

//     if (!logenc->get_active ()) {
//         filterPE.logenc         = falsePE.logenc;
//     }

//     if (!sh->get_active ()) {
//         filterPE.sh         = falsePE.sh;
//     }

//     if (!toneEqualizer->get_active()) {
//         filterPE.toneEqualizer = falsePE.toneEqualizer;
//     }

//     if (!epd->get_active ()) {
//         filterPE.epd        = falsePE.epd;
//     }

//     if (!fattal->get_active ()) {
//         filterPE.fattal     = falsePE.fattal;
//     }

//     if (!pcvignette->get_active ()) {
//         filterPE.pcvignette = falsePE.pcvignette;
//     }

//     if (!gradient->get_active ()) {
//         filterPE.gradient   = falsePE.gradient;
//     }

//     if (!labcurve->get_active ()) {
//         filterPE.labCurve   = falsePE.labCurve;
//     }

//     if (!sharpen->get_active ()) {
//         filterPE.sharpening      = falsePE.sharpening;
//     }

//     if (!sharpenmicro->get_active()) {
//         filterPE.sharpenMicro    = falsePE.sharpenMicro;
//     }

//     if (!impden->get_active ()) {
//         filterPE.impulseDenoise  = falsePE.impulseDenoise;
//     }

//     // if (!dirpyreq->get_active ()) {
//     //     filterPE.dirpyrequalizer = falsePE.dirpyrequalizer;
//     // }

//     if (!defringe->get_active ()) {
//         filterPE.defringe        = falsePE.defringe;
//     }

//     if (!dirpyrden->get_active ()) {
//         filterPE.denoise = falsePE.denoise;
//     }

//     if (!icm->get_active ()) {
//         filterPE.icm          = falsePE.icm;
//     }

//     if (!chmixer->get_active ()) {
//         filterPE.chmixer      = falsePE.chmixer;
//     }

//     if (!blackwhite->get_active ()) {
//         filterPE.blackwhite   = falsePE.blackwhite;
//     }

//     if (!filmSimulation->get_active ()) {
//         filterPE.filmSimulation  = falsePE.filmSimulation;
//     }

//     if (!softlight->get_active ()) {
//         filterPE.softlight = falsePE.softlight;
//     }

//     if (!dehaze->get_active ()) {
//         filterPE.dehaze = falsePE.dehaze;
//     }

//     if (!grain->get_active ()) {
//         filterPE.grain = falsePE.grain;
//     }

//     if (!rgbcurves->get_active ()) {
//         filterPE.rgbCurves    = falsePE.rgbCurves;
//     }

//     if (!distortion->get_active ()) {
//         filterPE.distortion   = falsePE.distortion;
//     }

//     if (!cacorr->get_active ()) {
//         filterPE.cacorrection = falsePE.cacorrection;
//     }

//     if (!vignetting->get_active ()) {
//         filterPE.vignetting   = falsePE.vignetting;
//     }

//     if (!lcp->get_active ()) {
//         filterPE.lensProf     = falsePE.lensProf;
//     }

//     if (!coarserot->get_active ()) {
//         filterPE.coarse      = falsePE.coarse;
//     }

//     if (!finerot->get_active ()) {
//         filterPE.rotate      = falsePE.rotate;
//     }

//     if (!crop->get_active ()) {
//         filterPE.crop        = falsePE.crop;
//     }

//     if (!resize->get_active ()) {
//         filterPE.resize      = falsePE.resize;
//     }

//     if (!prsharpening->get_active ()) {
//         filterPE.prsharpening      = falsePE.prsharpening;
//     }

//     if (!perspective->get_active ()) {
//         filterPE.perspective = falsePE.perspective;
//     }

//     if (!commonTrans->get_active ()) {
//         filterPE.commonTrans = falsePE.commonTrans;
//     }

//     if (!metadata->get_active()) {
//         filterPE.metadata = falsePE.metadata;
//     }

//     if (!exifch->get_active ()) {
//         filterPE.exif = falsePE.exif;
//     }

//     if (!iptc->get_active ()) {
//         filterPE.iptc = falsePE.iptc;
//     }

//     if (!colorcorrection->get_active()) {
//         filterPE.colorcorrection = falsePE.colorcorrection;
//     }
    
//     if (!smoothing->get_active()) {
//         filterPE.smoothing = falsePE.smoothing;
//     }
    
//     if (!cbdl->get_active()) {
//         filterPE.dirpyrequalizer = falsePE.dirpyrequalizer;
//     }    
    
//     if (!raw_method->get_active ()) {
//         filterPE.raw.bayersensor.method   = falsePE.raw.bayersensor.method;
//         filterPE.raw.bayersensor.dualDemosaicAutoContrast = falsePE.raw.bayersensor.dualDemosaicAutoContrast;
//         filterPE.raw.bayersensor.dualDemosaicContrast = falsePE.raw.bayersensor.dualDemosaicContrast;
//         filterPE.raw.xtranssensor.method  = falsePE.raw.xtranssensor.method;
//         filterPE.raw.xtranssensor.dualDemosaicAutoContrast = falsePE.raw.xtranssensor.dualDemosaicAutoContrast;
//         filterPE.raw.xtranssensor.dualDemosaicContrast = falsePE.raw.xtranssensor.dualDemosaicContrast;
//     }

//     if (!raw_border->get_active ()) {
//         filterPE.raw.bayersensor.border = falsePE.raw.bayersensor.border;
//         filterPE.raw.xtranssensor.border = falsePE.raw.xtranssensor.border;
//     }

//     if (!raw_imagenum->get_active ()) {
//         filterPE.raw.bayersensor.imageNum = falsePE.raw.bayersensor.imageNum;
//     }

//     if (!raw_ccSteps->get_active ()) {
//         filterPE.raw.bayersensor.ccSteps  = falsePE.raw.bayersensor.ccSteps;
//         filterPE.raw.xtranssensor.ccSteps = falsePE.raw.xtranssensor.ccSteps;
//     }

//     if (!raw_dcb_iterations->get_active ()) {
//         filterPE.raw.bayersensor.dcbIterations   = falsePE.raw.bayersensor.dcbIterations;
//     }

//     if (!raw_dcb_enhance->get_active ()) {
//         filterPE.raw.bayersensor.dcbEnhance      = falsePE.raw.bayersensor.dcbEnhance;
//     }

//     //if (!raw_all_enhance->get_active ())     filterPE.raw.bayersensor.allEnhance      = falsePE.raw.bayersensor.allEnhance;
//     if (!raw_lmmse_iterations->get_active ()) {
//         filterPE.raw.bayersensor.lmmseIterations = falsePE.raw.bayersensor.lmmseIterations;
//     }

//     if (!raw_black->get_active ()) {
//         filterPE.raw.bayersensor.exBlack0        = falsePE.raw.bayersensor.exBlack0;
//         filterPE.raw.bayersensor.exBlack1        = falsePE.raw.bayersensor.exBlack1;
//         filterPE.raw.bayersensor.exBlack2        = falsePE.raw.bayersensor.exBlack2;
//         filterPE.raw.bayersensor.exBlack3        = falsePE.raw.bayersensor.exBlack3;
//         filterPE.raw.bayersensor.exTwoGreen      = falsePE.raw.bayersensor.exTwoGreen;
//         filterPE.raw.xtranssensor.exBlackRed     = falsePE.raw.xtranssensor.exBlackRed;
//         filterPE.raw.xtranssensor.exBlackGreen   = falsePE.raw.xtranssensor.exBlackGreen;
//         filterPE.raw.xtranssensor.exBlackBlue    = falsePE.raw.xtranssensor.exBlackBlue;
//     }

//     if (!raw_pixelshift->get_active ()) {
//         filterPE.raw.bayersensor.pixelShiftBlur                   = falsePE.raw.bayersensor.pixelShiftBlur;
//         filterPE.raw.bayersensor.pixelShiftEperIso                = falsePE.raw.bayersensor.pixelShiftEperIso;
//         filterPE.raw.bayersensor.pixelShiftEqualBright            = falsePE.raw.bayersensor.pixelShiftEqualBright;
//         filterPE.raw.bayersensor.pixelShiftEqualBrightChannel     = falsePE.raw.bayersensor.pixelShiftEqualBrightChannel;
//         filterPE.raw.bayersensor.pixelShiftGreen                  = falsePE.raw.bayersensor.pixelShiftGreen;
//         filterPE.raw.bayersensor.pixelShiftHoleFill               = falsePE.raw.bayersensor.pixelShiftHoleFill;
//         filterPE.raw.bayersensor.pixelShiftDemosaicMethod         = falsePE.raw.bayersensor.pixelShiftDemosaicMethod;
//         filterPE.raw.bayersensor.pixelShiftMedian                 = falsePE.raw.bayersensor.pixelShiftMedian;
//         filterPE.raw.bayersensor.pixelShiftMotionCorrectionMethod = falsePE.raw.bayersensor.pixelShiftMotionCorrectionMethod;
//         filterPE.raw.bayersensor.pixelShiftNonGreenCross          = falsePE.raw.bayersensor.pixelShiftNonGreenCross;
//         filterPE.raw.bayersensor.pixelShiftSigma                  = falsePE.raw.bayersensor.pixelShiftSigma;
//         filterPE.raw.bayersensor.pixelShiftSmooth                 = falsePE.raw.bayersensor.pixelShiftSmooth;
//         filterPE.raw.bayersensor.pixelShiftShowMotion             = falsePE.raw.bayersensor.pixelShiftShowMotion;
//         filterPE.raw.bayersensor.pixelShiftShowMotionMaskOnly     = falsePE.raw.bayersensor.pixelShiftShowMotionMaskOnly;
//     }

//     if (!raw_linenoise->get_active ()) {
//         filterPE.raw.bayersensor.linenoise       = falsePE.raw.bayersensor.linenoise;
//         filterPE.raw.bayersensor.linenoiseDirection = falsePE.raw.bayersensor.linenoiseDirection;
//     }

//     if (!raw_greenthresh->get_active ()) {
//         filterPE.raw.bayersensor.greenEq         = falsePE.raw.bayersensor.greenEq;
//     }

//     if (!raw_expos->get_active ()) {
//         filterPE.raw.exPos              = falsePE.raw.exPos;
//     }

//     if (!raw_ca_autocorrect->get_active ()) {
//         filterPE.raw.ca_autocorrect       = falsePE.raw.ca_autocorrect;
//         filterPE.raw.caautoiterations       = falsePE.raw.caautoiterations;
//     }

//     if (!raw_caredblue->get_active ()) {
//         filterPE.raw.cared              = falsePE.raw.cared;
//         filterPE.raw.cablue             = falsePE.raw.cablue;
//     }

//     if (!raw_ca_avoid_colourshift->get_active ()) {
//         filterPE.raw.ca_avoidcolourshift = falsePE.raw.ca_avoidcolourshift;
//     }

//     if (!raw_hotpix_filt->get_active ())     {
//         filterPE.raw.hotPixelFilter     = falsePE.raw.hotPixelFilter;
//     }

//     if (!raw_deadpix_filt->get_active ())    {
//         filterPE.raw.deadPixelFilter    = falsePE.raw.deadPixelFilter;
//     }

//     if (!raw_deadpix_filt->get_active () && !raw_hotpix_filt->get_active ()) {
//         filterPE.raw.hotdeadpix_thresh = falsePE.raw.hotdeadpix_thresh;
//     }

//     if (!raw_pdaf_lines_filter->get_active ())    {
//         filterPE.raw.bayersensor.pdafLinesFilter = falsePE.raw.bayersensor.pdafLinesFilter;
//     }

//     if (!df_file->get_active ()) {
//         filterPE.raw.darkFrame          = falsePE.raw.darkFrame;
//     }

//     if (!df_AutoSelect->get_active ()) {
//         filterPE.raw.df_autoselect             = falsePE.raw.df_autoselect;
//     }

//     if (!ff_file->get_active ()) {
//         filterPE.raw.ff_file            = falsePE.raw.ff_file;
//     }

//     if (!ff_AutoSelect->get_active ()) {
//         filterPE.raw.ff_AutoSelect      = falsePE.raw.ff_AutoSelect;
//     }

//     if (!ff_BlurRadius->get_active ()) {
//         filterPE.raw.ff_BlurRadius      = falsePE.raw.ff_BlurRadius;
//     }

//     if (!ff_BlurType->get_active ()) {
//         filterPE.raw.ff_BlurType        = falsePE.raw.ff_BlurType;
//     }

//     if (!ff_ClipControl->get_active ()) {
//         filterPE.raw.ff_clipControl     = falsePE.raw.ff_clipControl;
//         filterPE.raw.ff_AutoClipControl = falsePE.raw.ff_AutoClipControl;
//     }

//     if (dstPE) {
//         *dstPE = filterPE;
//     }

//     // Apply the filter!
//     filterPE.combine(*dstPP, *srcPP, true);
// }


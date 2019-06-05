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
#include "toneequalizer.h"
#include "eventmapper.h"

using namespace rtengine;
using namespace rtengine::procparams;


ToneEqualizer::ToneEqualizer(): FoldableToolPanel(this, "toneequalizer", M("TP_TONE_EQUALIZER_LABEL"), false, true)
{
    auto m = ProcEventMapper::getInstance();
    EvEnabled = m->newEvent(RGBCURVE, "HISTORY_MSG_TONE_EQUALIZER_ENABLED");
    EvBands = m->newEvent(RGBCURVE, "HISTORY_MSG_TONE_EQUALIZER_BANDS");
    EvDetail = m->newEvent(RGBCURVE, "HISTORY_MSG_TONE_EQUALIZER_DETAIL");

    for (size_t i = 0; i < bands.size(); ++i) {
        bands[i] = Gtk::manage(new Adjuster(M("TP_TONE_EQUALIZER_BAND_" + std::to_string(i)), -100, 100, 1, 0));
        bands[i]->setAdjusterListener(this);
        pack_start(*bands[i]);
    }
    pack_start(*Gtk::manage(new Gtk::HSeparator()));
    detail = Gtk::manage(new Adjuster(M("TP_TONE_EQUALIZER_DETAIL"), -5, 5, 1, 0));
    detail->setAdjusterListener(this);
    pack_start(*detail);
    
    show_all_children ();
}


void ToneEqualizer::read(const ProcParams *pp)
{
    disableListener();

    setEnabled(pp->toneEqualizer.enabled);

    for (size_t i = 0; i < bands.size(); ++i) {
        bands[i]->setValue(pp->toneEqualizer.bands[i]);
    }
    detail->setValue(pp->toneEqualizer.detail);
    
    enableListener();
}


void ToneEqualizer::write(ProcParams *pp)
{
    for (size_t i = 0; i < bands.size(); ++i) {
        pp->toneEqualizer.bands[i] = bands[i]->getValue();
    }
    pp->toneEqualizer.enabled = getEnabled();
    pp->toneEqualizer.detail = detail->getValue();
}


void ToneEqualizer::setDefaults(const ProcParams *defParams)
{
    for (size_t i = 0; i < bands.size(); ++i) {
        bands[i]->setDefault(defParams->toneEqualizer.bands[i]);
    }
    detail->setDefault(defParams->toneEqualizer.detail);
}


void ToneEqualizer::adjusterChanged(Adjuster *a, double newval)
{
    if (listener && getEnabled()) {
        if (a == detail) {
            listener->panelChanged(EvDetail, Glib::ustring::format(a->getValue()));
        } else {
            Glib::ustring s;
            for (size_t i = 0; i < bands.size(); ++i) {
                s += Glib::ustring::format((int)bands[i]->getValue()) + " ";
            }
            listener->panelChanged(EvBands, s);
        }
    }
}


void ToneEqualizer::adjusterAutoToggled(Adjuster *a, bool newval)
{
}


void ToneEqualizer::enabledChanged()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged(EvEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged(EvEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvEnabled, M("GENERAL_DISABLED"));
        }
    }
}


void ToneEqualizer::trimValues(rtengine::procparams::ProcParams *pp)
{
    for (size_t i = 0; i < bands.size(); ++i) {
        bands[i]->trimValue(pp->toneEqualizer.bands[i]);
    }
    detail->trimValue(pp->toneEqualizer.detail);
}

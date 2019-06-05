/* -*- C++ -*-
 *  
 *  This file is part of RawTherapee.
 */
#ifndef _PCVIGNETTE_H_
#define _PCVIGNETTE_H_

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"

class PCVignette : public ToolParamBlock, public AdjusterListener, public FoldableToolPanel
{

protected:
    Adjuster* strength;
    Adjuster* feather;
    Adjuster* roundness;

public:

    PCVignette ();

    void read(const rtengine::procparams::ProcParams* pp) override;
    void write(rtengine::procparams::ProcParams* pp) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams) override;
    void adjusterChanged (Adjuster* a, double newval) override;
    void adjusterAutoToggled(Adjuster* a, bool newval) override;
    void enabledChanged  () override;
    void trimValues          (rtengine::procparams::ProcParams* pp) override;
};

#endif

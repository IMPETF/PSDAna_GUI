#include "TApplication.h"
#include "GuiFrame.h"
#include "TGClient.h"

int main(int argc, char *argv[])
{

    TApplication apponline("PSDAna_GUI", &argc, argv);

    GuiFrame *mygui = new GuiFrame(gClient->GetRoot(), 10, 10);

    apponline.Run();

    return 0;
}
